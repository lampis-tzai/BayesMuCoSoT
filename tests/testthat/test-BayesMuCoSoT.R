test_that("random data test", {
  all_data = iris[iris$Species=='setosa',]
  questioned_data = all_data[1:(nrow(all_data)/2),]
  known_data = all_data[(nrow(all_data)/2+1):nrow(all_data),]
  background_data = iris[iris$Species!='setosa',]
  background_data_id = 'Species'
  y = names(questioned_data)[1:4]
  expect_equal(BayesMuCoSoT_fit(y,x=NA,questioned_data,known_data,
                                background_data,background_data_id), 20,
                                tolerance=2)
})
