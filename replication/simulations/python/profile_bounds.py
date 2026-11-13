import heapq
import itertools
import numpy as np
from scipy.optimize import minimize


class Profile:
    def __init__(self, mu, coefficients, upper, dimensions):
        self.mu = [np.asarray(x, dtype=float) for x in mu]
        self.coef = np.asarray(coefficients, dtype=float)
        self.upper = np.asarray(upper, dtype=float)
        self.dimensions = np.asarray(dimensions, dtype=int)
        self.order = len(self.mu)
        self.active = np.flatnonzero(self.upper > 0)
        self.masks = (
            np.arange(2**self.order)[:, None] >> np.arange(self.order) & 1
        ).astype(bool)
        self.corners = self.masks.astype(float)
        if (
            self.coef.shape != (2**self.order,)
            or self.upper.shape != (self.order,)
            or self.dimensions.shape != (self.order,)
            or (not np.all(np.isfinite(self.coef)))
            or (not np.all(np.isfinite(self.upper)))
            or np.any(self.upper < 0)
            or np.any(self.dimensions < 1)
        ):
            raise ValueError("invalid fixed-profile components")
        for m, eigen in enumerate(self.mu):
            if (
                not np.all(np.isfinite(eigen))
                or (self.upper[m] > 0 and len(eigen) != self.dimensions[m])
                or np.any(1 + self.upper[m] * eigen <= 0)
            ):
                raise ValueError("invalid spectral domain")
        if np.min(self.traces(self.corners * self.upper)) <= 0:
            raise ValueError("trace is not positive throughout the box")

    def traces(self, points):
        points = np.atleast_2d(np.asarray(points, dtype=float))
        terms = np.prod(
            np.where(self.masks[None, :, :], points[:, None, :], 1.0), axis=2
        )
        return terms @ self.coef

    def trace(self, point):
        return float(self.traces(point)[0])

    def value(self, point):
        point = np.asarray(point, dtype=float)
        determinant = sum(
            (
                np.log1p(point[m] * eigen).sum() / self.dimensions[m]
                for m, eigen in enumerate(self.mu)
            )
        )
        return float(determinant - np.log(self.trace(point) / self.coef[0]))

    def gradient(self, point):
        point = np.asarray(point, dtype=float)
        trace = self.trace(point)
        gradient = np.zeros(self.order)
        for m in self.active:
            select = self.masks[:, m]
            others = np.arange(self.order) != m
            terms = np.prod(
                np.where(self.masks[select][:, others], point[others], 1.0), axis=1
            )
            derivative = terms @ self.coef[select]
            eigen = self.mu[m]
            gradient[m] = (
                np.sum(eigen / (1 + point[m] * eigen)) / self.dimensions[m]
                - derivative / trace
            )
        return gradient

    def box_upper(self, lower, upper):
        lower, upper = (np.asarray(lower), np.asarray(upper))
        centre = (lower + upper) / 2
        vertices = lower + self.corners * (upper - lower)
        ratios = np.ones(len(vertices))
        for m, eigen in enumerate(self.mu):
            if len(eigen):
                terms = (1 + vertices[:, m, None] * eigen) / (1 + centre[m] * eigen)
                ratios *= terms.mean(axis=1)
        ratios *= self.trace(centre) / self.traces(vertices)
        maximum = float(np.max(ratios))
        value = self.value(centre)
        padding = 256 * np.finfo(float).eps * (1 + abs(value) + abs(maximum))
        return value + max(0.0, maximum - 1.0) + padding

    def maximise(self, tolerance, max_boxes=20000, seed=None):
        if tolerance <= 0 or max_boxes < 1:
            raise ValueError("positive tolerance and box budget are required")
        if not len(self.active):
            return dict(
                lower=0.0,
                upper=0.0,
                point=[0.0] * self.order,
                resolved=True,
                boxes=0,
                starts=0,
                successes=0,
            )
        best_point = np.zeros(self.order)
        best = self.value(best_point)
        starts = list(itertools.product(np.linspace(0, 1, 5), repeat=len(self.active)))
        if seed is not None:
            seed = np.asarray(seed, dtype=float)
            if (
                seed.shape != self.upper.shape
                or np.any(seed < 0)
                or np.any(seed > self.upper)
            ):
                raise ValueError("infeasible replayed solution")
            best_point, best = (seed.copy(), self.value(seed))
            starts.append(tuple(seed[self.active] / self.upper[self.active]))

        def expand(z):
            point = np.zeros(self.order)
            point[self.active] = z * self.upper[self.active]
            return point

        def objective(z):
            point = expand(z)
            return (
                -self.value(point),
                -self.gradient(point)[self.active] * self.upper[self.active],
            )

        successes = 0
        for start in starts:
            fit = minimize(
                objective,
                start,
                method="L-BFGS-B",
                jac=True,
                bounds=[(0, 1)] * len(self.active),
                options=dict(ftol=1e-15, gtol=1e-10, maxiter=1000, maxls=40),
            )
            successes += int(fit.success)
            point = expand(np.clip(fit.x, 0, 1))
            value = self.value(point)
            if np.isfinite(value) and value > best:
                best_point, best = (point, value)
        lower, upper = (np.zeros(self.order), self.upper.copy())
        queue = [(-self.box_upper(lower, upper), 0, lower, upper)]
        boxes = 1
        serial = 0
        while queue and -queue[0][0] - best > tolerance and (boxes + 2 <= max_boxes):
            _, _, lower, upper = heapq.heappop(queue)
            widths = (upper - lower) / np.where(self.upper > 0, self.upper, 1)
            axis = int(np.argmax(widths))
            middle = (lower[axis] + upper[axis]) / 2
            if middle in (lower[axis], upper[axis]):
                raise ArithmeticError("subdivision reached floating-point resolution")
            left_upper, right_lower = (upper.copy(), lower.copy())
            left_upper[axis] = right_lower[axis] = middle
            for lo, hi in ((lower, left_upper), (right_lower, upper)):
                centre = (lo + hi) / 2
                value = self.value(centre)
                if value > best:
                    best_point, best = (centre, value)
                bound = self.box_upper(lo, hi)
                serial += 1
                boxes += 1
                if bound >= best:
                    heapq.heappush(queue, (-bound, serial, lo, hi))
        bound = max(best, -queue[0][0]) if queue else best
        return dict(
            lower=float(best),
            upper=float(bound),
            point=best_point.tolist(),
            resolved=bool(bound - best <= tolerance),
            boxes=boxes,
            starts=len(starts),
            successes=successes,
        )
