# Copy to config.sh and edit on MSI before submitting setup.
MSI_ACCOUNT="REPLACE_WITH_YOUR_MSI_ACCOUNT"
# Find actual compatible versions with `module spider R` and `module spider gsl`.
# Include the compiler module if your R/GSL modules require it. R >= 4.4 is needed.
MSI_MODULES=("REPLACE_WITH_R_MODULE" "REPLACE_WITH_GSL_MODULE")
MSI_PARTITION="msismall"
MSI_TIME="24:00:00"
MSI_MEMORY="16G"
MSI_CONCURRENCY=4
# Prefer a persistent project directory with sufficient quota (reserve >= 30 GB).
# ROOT is set by the launcher to the extracted bundle's absolute directory.
MSI_RESULTS="$ROOT/runs"
MSI_LIBRARY="$ROOT/library-linux"
