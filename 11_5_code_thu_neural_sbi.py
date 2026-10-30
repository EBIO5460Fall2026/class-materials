# Neural Simulation Based Inference (SBI)
# Demonstrate how a neural network can predict parameters from data
# Train a neural network

# To create this computing environment in conda:
# conda create --name py-tensorflow-pm4e26
# conda activate py-tensorflow-pm4e26
# conda install python=3.12 -c conda-forge
# pip install tensorflow numpy pandas matplotlib

import numpy as np
import matplotlib.pyplot as plt
import keras                      #neural network library, keras 3.15.1
import tensorflow as tf           #neural network library, tensorflow 2.21.0

# Warnings and CUDA errors from keras and tensorflow are harmless (this merely
# means you don't have an NVIDIA GPU).

# First load our large training set of time series labeled with parameter values
# sampled from the prior distribution. These were made in the previous script.

c_sample = np.load("c_sample.npy")
N_sim_grid = np.load("N_sim_grid.npy")


# Neural net to learn c from trajectories

keras.utils.set_random_seed(6590)

x_train = N_sim_grid       # shape (50000, 50)
y_train = c_sample * 1000  # c values (scaled to about 1)

# Neural network
# A standard regression network
# CUDA error is harmless
model = keras.Sequential([
    keras.layers.Input(shape=(x_train.shape[1],)),
    keras.layers.Dense(64, activation="relu"),
    keras.layers.Dense(64, activation="relu"),
    keras.layers.Dense(1)
])

model.summary()

model.compile(
    optimizer="rmsprop",
    loss="mse"
)

history = model.fit(
    x_train,
    y_train,
    epochs=30,
    batch_size=32,
    validation_split=0.2
)

# Training is mostly done after 30 epochs:

plt.plot(history.history["loss"], label="training")
plt.plot(history.history["val_loss"], label="validation")
plt.xlabel("Epoch")
plt.ylabel("MSE")
plt.ylim(bottom=0.1, top=0.2)
plt.legend()

# Now we have a trained network, can it predict c for our dataset?

N_obs_grid = np.load("N_obs_grid.npy")
x_test = N_obs_grid.reshape(1, 50) #correct shape to enter network
c_pred = model.predict(x_test) / 1000

# It predicts rather well
print(f"Predicted c: {c_pred}, True c: 0.001")

# The above is predicting c from y, i.e.
#     predict c | y
# which is just a point estimate. We want it to predict the posterior
# probabilities for different values of c, i.e.
#     predict p(c | y)
# We do that next.


# Neural net to learn p(c | y) from trajectories

keras.utils.set_random_seed(76801)

x_train = N_sim_grid       # shape (50000, 50), unchanged from above

# First bin c into say 30 equal classes. These will correspond to histogram bars
# of the posterior.
# Classes: 0, 1, ..., 29

n_bins = 30
bin_edges = np.linspace(c_sample.min(), c_sample.max(), n_bins + 1)
y_train = np.digitize(c_sample, bin_edges[1:-1])

# Neural network
# A standard classification network with a softmax final layer
model = keras.Sequential([
    keras.layers.Input(shape=(x_train.shape[1],)),
    keras.layers.Dense(64, activation="relu"),
    keras.layers.Dense(64, activation="relu"),
    keras.layers.Dense(n_bins, activation="softmax")
])

model.summary()

# Cross entropy loss
model.compile(
    optimizer="rmsprop",
    loss="sparse_categorical_crossentropy",
    metrics=["accuracy"]
)

history = model.fit(
    x_train,
    y_train,
    epochs=50,
    batch_size=32,
    validation_split=0.2
)

# Training is mostly done after 35 epochs and then perhaps gets worse:

plt.figure()
plt.plot(history.history["loss"], label="training")
plt.plot(history.history["val_loss"], label="validation loss")
plt.xlabel("Epoch")
plt.ylabel("Crossentropy loss")
plt.legend()

# Accuracy is lowish and not improving:

plt.figure()
plt.plot(history.history["val_accuracy"], label="validation accuracy")
plt.xlabel("Epoch")
plt.ylabel("Accuracy")
plt.ylim(bottom=0, top=1)
plt.legend()

# but we'll see next this is nevertheless a good picture of the posterior.

# What is the predicted p(c|y) for our dataset?

c_post_pred = model.predict(x_test)
c_post_pred.sum() #sum to 1
c_post_pred.shape
c_post_pred = c_post_pred[0]

c_true = 0.001
plt.figure()
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
plt.bar(bin_centers, c_post_pred, width=np.diff(bin_edges), align="center")
plt.axvline(c_true, color="orange", label="true value")
plt.xlabel("c")
plt.ylabel("Posterior probability")
plt.legend()

# This looks very similar to the result from rejection ABC