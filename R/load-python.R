py_optimizer <- reticulate::import_from_path("py_optimizer", file.path("inst", "python"))

for (obj in names(py_optimizer)) {
  assign(obj, NULL)
}
# Clean up
rm(py_optimizer)

# Now all those names are in the namespace, and ready to be replaced on load
.onLoad <- function(libname, pkgname) {

    py_optimizer <- reticulate::import_from_path(
                                    "py_optimizer",
                                    system.file("python", package = packageName()))

  ## assignInMyNamespace(...) is meant for namespace manipulation
  for (obj in names(py_optimizer)) {
    assignInMyNamespace(obj, py_optimizer[[obj]])
  }
}
