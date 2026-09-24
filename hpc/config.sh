# This configuration describes Anvil resources; scientific controls live in R.
IMR_ACCOUNT=cis260924
IMR_PARTITION=shared
IMR_CONCURRENCY=16
IMR_MODULES=(gcc/11.2.0 r/4.4.1 gsl/2.4 python/3.9.5)
IMR_DEPENDENCIES="$IMR_ROOT/runtime/dependencies"
IMR_RUNS="$IMR_ROOT/runs"
