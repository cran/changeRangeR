#' @name complementarity
#' @title Compare raster values within and outside of a mask
#' @param ras1 `raster` object of categorical values (e.g., a species richness map)
#' @param ras1mask `raster` object representing areas of interest (e.g., protected areas)
#' @return A list of two objects. The first is the proportion of the range within the mask, and the second is the proportion of unique values masked.
#' @examples
#' # create raster
#' ras1 <- raster::raster(nrows=108, ncols=108, xmn=-50, xmx=50)
#' raster::values(ras1)<- runif(n = (108*108))
#' ras1[ras1 < 0.5] <- NA
#' ras1[!is.na(ras1)] <- 1
#' # create ras1mask
#' ras1mask <- raster::raster(nrows=108, ncols=108, xmn=-50, xmx=50)
#' raster::values(ras1mask)<- runif(n = (108*108))
#' ras1mask[ras1mask < 0.15] <- NA
#' ras1mask[!is.na(ras1mask)] <- 1
#' # complementarity
#' complementarity(ras1, ras1mask)
#' @export
#'

complementarity <- function(ras1, ras1mask){
  #library(raster)
  ## calculate percentage of values protected
  rasMask <- raster::mask(ras1, ras1mask)
  prop <- raster::cellStats(rasMask, stat = "sum", na.rm = T) / raster::cellStats(ras1, stat = "sum", na.rm = T) * 100
  ## Get percentage of values that fall in mask
  SRdf <- as.data.frame(table(raster::values(ras1)))
  df <- as.data.frame(table(raster::values(rasMask)))
  # merge datafames
  dfMerged <- merge(SRdf, df, by = "Var1", all = T)
  dfMerged[is.na(dfMerged)] <- 0
  dfMerged$percent <- dfMerged$Freq.y / dfMerged$Freq.x * 100
  # Combine outputs
  out <- list(Percent_of_Total = prop, Percent_unique_values = dfMerged)
  return(out)
}

