test_that("same source test", {
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

test_that("different source test", {
  set.seed(10)
  questioned_data = iris[iris$Species=='versicolor',]
  split_ind <- sample(seq_len(nrow(questioned_data)), size = floor(0.5 * nrow(questioned_data)))
  background_data1 = questioned_data[-split_ind,]
  questioned_data = questioned_data[split_ind,]


  known_data = iris[iris$Species=='virginica',]
  split_ind <- sample(seq_len(nrow(known_data)), size = floor(0.5 * nrow(known_data)))
  background_data2 = known_data[-split_ind,]
  known_data = known_data[split_ind,]

  background_data = rbind(iris[iris$Species=='setosa',],
                          background_data1,
                          background_data2)

  background_data_id = 'Species'
  y = names(questioned_data)[1:4]
  expect_equal(BayesMuCoSoT_fit(y,x=NA,questioned_data,known_data,
                                background_data,background_data_id), -30,
               tolerance=2)
})
