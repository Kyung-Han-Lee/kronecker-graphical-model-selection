N_SWEEP <- 2L

SUPPORT_TOL <- 1e-6

SQRT_METHOD <- "none"

GAMMA_GRID <- c(0, 0.25, 0.5, 1)

JOINT_U_CAP <- 10

JOINT_BOUNDARY_MARGIN <- 1e-8

ETA_LO <- 1e-4

ETA_HI <- 100

N_ETA_DEFAULT <- 40L

ARM_SPECS <- list(
    control = list(
        subdir = "combat", label = "TR=3s, p_T=187", p_time = 187L,
        requires_manifest = FALSE
    ),
    resampled = list(
        subdir = "combat_tr5", label = "TR=5s, p_T=113", p_time = 113L,
        requires_manifest = TRUE
    ),
    projected = list(

        subdir = "combat_pb102", label = "통과대역 사영, p=102", p_time = 102L,
        requires_manifest = TRUE
    )
)

GAMMA_REPORT <- 0.5

BOOT_PAIRS <- list(
    posterior_cingulate = c(35L, 36L),
    precuneus = c(67L, 68L)
)

GROUP_SEED_OFFSET <- c(CN = 1L, MCI = 2L, Dementia = 3L)

make_eta_grid <- function(n_eta) {
    stopifnot(length(n_eta) == 1L, is.finite(n_eta), n_eta >= 3L)
    exp(seq(log(ETA_LO), log(ETA_HI), length.out = n_eta))
}
