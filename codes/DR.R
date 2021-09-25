DR <- function (para = list(train_data,
                            test_data,
                            denoising_rate = 0.2, 
                            lambda = 1e-5,
                            activation_hidden = 'sigmoid',
                            activation_output = 'sigmoid',
                            batch_size = 256,
                            learning_rate = 0.001,
                            epochs = 40,
                            encoded_1 = 1024,
                            encoded_2 = 128)){
  
  library("keras")
  library("caret")
  
  train_data <- para[['train_data']]
  test_data <- para[['test_data']]
  denoising_rate <- para[['denoising_rate']]
  lambda <- para[['lambda']]
  activation_hidden <- para[['activation_hidden']]
  activation_output <- para[['activation_output']]
  batch_size <- para[['batch_size']]
  learning_rate <- para[['learning_rate']]
  epochs <- para[['epochs']]
  encoded_1 <- para[['encoded_1']]
  encoded_2 <- para[['encoded_2']]
  
  
  microtrain_layer1 <- train_data
  microtest_layer1 <- test_data
  final_proj_dim <- ncol(train_data)
  
  keras::k_clear_session()
  
  tns = encoder_input = keras::layer_input(shape = final_proj_dim)
  tns = encoder_denoise = keras::layer_dropout(object = encoder_input, rate = denoising_rate)
  
  tns = keras::layer_dense(tns, units = encoded_1, activation = activation_hidden, activity_regularizer = keras::regularizer_l1(lambda))
  tns = keras::layer_activation_leaky_relu(tns)
  encoder = keras::keras_model(inputs = encoder_input, outputs = tns)
  
  tns = decoder_input = keras::layer_input(shape = encoded_1)
  tns = keras::layer_dense(tns, units = final_proj_dim, activation = activation_output, activity_regularizer = keras::regularizer_l1(lambda))
  decoder = keras::keras_model(inputs = decoder_input, outputs = tns)
  
  tns = ae_input = keras::layer_input(final_proj_dim)
  tns = decoder(encoder(tns))
  autoencoder = keras::keras_model(inputs = ae_input, outputs = tns)
  keras::compile(autoencoder, optimizer = keras::optimizer_adam(lr = learning_rate), loss = 'mean_squared_error')
  
  keras::fit(autoencoder, microtrain_layer1, microtrain_layer1, batch_size = batch_size, epochs = epochs, verbose = 0)
  
  microtrain_layer2 <- predict(encoder, microtrain_layer1, batch_size = batch_size)
  microtest_layer2 <- predict(encoder, microtest_layer1, batch_size = batch_size)
  
  tns = encoder_input = keras::layer_input(shape = encoded_1)
  tns = encoder_denoise = keras::layer_dropout(object = encoder_input, rate = denoising_rate)
  tns = keras::layer_dense(tns, units = encoded_2, activation = activation_hidden, activity_regularizer = keras::regularizer_l1(lambda))
  tns = keras::layer_activation_leaky_relu(tns)
  encoder = keras::keras_model(inputs = encoder_input, outputs = tns)
  
  tns = decoder_input = keras::layer_input(shape = encoded_2)
  tns = keras::layer_dense(tns, units = encoded_1, activation = activation_output, activity_regularizer = keras::regularizer_l1(lambda))
  decoder = keras::keras_model(inputs = decoder_input, outputs = tns)
  
  tns = ae_input= keras::layer_input(encoded_1)
  tns = decoder(encoder(tns))
  autoencoder = keras::keras_model(inputs = ae_input, outputs = tns)
  keras::compile(autoencoder, optimizer = keras::optimizer_adam(lr = learning_rate), loss = 'mean_squared_error')
  
  keras::fit(autoencoder, microtrain_layer2, microtrain_layer2, batch_size = batch_size, epochs = epochs, verbose = 0)
  
  microtrain_data <- predict(encoder, microtrain_layer2, batch_size = batch_size)
  microtest_data <- predict(encoder, microtest_layer2, batch_size = batch_size)
  return(list(microtrain_data,microtest_data))
}