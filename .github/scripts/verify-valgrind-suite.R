# Both parent and children run under Memcheck. An instrumented parent compared
# with ordinary children is not a same-execution-environment numerical check.
scripts <- file.path(Sys.getenv("GITHUB_WORKSPACE"), ".github", "scripts")
source(file.path(scripts, "configure-valgrind-workers.R"))
source(file.path(scripts, "verify-valgrind-logs.R"))
source(file.path(scripts, "verify-tests.R"))
verify_valgrind_worker_logs()
