# create continuous raster
p <- raster::raster(nrows=108, ncols=108, xmn=-50, xmx=50)
raster::values(p)<- runif(n = (108*108))
# create occurrences
xy <- dismo::randomPoints(p, 4)
# create original convex hull
ch.orig <- mcp(xy)
# set threshold
thr <- 0.5
# mcpSDM
out <- mcpSDM(p, xy, ch.orig, thr)


## TEST
test_that("output type checks", {
  expect_type(out, "list")
  expect_is(out$jsi, "numeric")
  expect_is(out$thr, "numeric")
  expect_is(out$ov.pts, "numeric")
  expect_is(out$best.fit, "SpatialPolygons")
  expect_is(out$best.fit.ind, "integer")
})
