#' @export prepareBIOM


prepareBIOM <- function(
  asvtab,
  write = TRUE,
  filename = "bacteria.biom"
){
  result <- biomformat::make_biom(data = t(asvtab))
  if(write){
    biomformat::write_biom(result, filename)
  }

  return(result)
}
