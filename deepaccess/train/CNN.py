import keras
import numpy as np
from keras.models import Sequential
from keras.layers import (
    Conv1D,
    Dense,
    Dropout,
    LSTM,
    GlobalMaxPooling1D,
    MaxPooling1D,
    Flatten,
)
from keras import optimizers
import pkg_resources

class CNN:
    # CNN code for motif convolutional layer
    # taken from https://github.com/uci-cbcl/DanQ
    def __init__(
        self,
        model_layers,
        seq_shape,
        out_shape=2,
        conv_filter_number=100,
        conv_filter_size=20,
    ):
        model = Sequential()
        for i, layer in enumerate(model_layers):
            if layer == "conv":
                if i == 0:
                    HOMER_motifs = list(
                        np.load(
                            pkg_resources.resource_stream(
                                __name__,
                                "homer_matrix.npy"
                            ),
                            allow_pickle=True))
                    conv_layer = Conv1D(
                        input_shape=seq_shape,
                        filters=len(HOMER_motifs) * 2,
                        kernel_size=35,
                        padding="valid",
                        activation="relu",
                        strides=1,
                    )
                    model.add(conv_layer)
                    conv_weights = conv_layer.get_weights()

                    reverse_motifs = [
                        HOMER_motifs[j][::-1, ::-1]
                        for j in range(len(HOMER_motifs))
                    ]
                    HOMER_motifs = HOMER_motifs + reverse_motifs

                    for j in range(len(HOMER_motifs)):
                        m = HOMER_motifs[j][::-1, :]
                        w = m.shape[0]
                        conv_weights[0][:, :, j] = 0
                        start = np.random.randint(low=3, high=35 - w - 3 + 1)
                        conv_weights[0][start:(start + w), :, j] = m - 0.25
                        conv_weights[1][j] = np.random.uniform(low=-1.0,
                                                               high=0.0)

                    conv_layer.set_weights(conv_weights)
                    conv_layer.trainable = False
                else:
                    model.add(
                        Conv1D(
                            conv_filter_number,
                            conv_filter_size,
                            activation="relu",
                            padding="same",
                        )
                    )
            if layer == "globalpool":
                model.add(GlobalMaxPooling1D())
            if layer == "maxpool":
                model.add(MaxPooling1D(3, 1))
            if layer == "LSTM":
                model.add(LSTM(16, return_sequences=True))
            if layer == "dense":
                model.add(Dropout(0.1))
                model.add(Dense(128, activation="relu"))
        if "dense" not in model_layers and "globalpool" not in model_layers:
            model.add(Flatten())

        model.add(Dropout(0.1))
        model.add(Dense(out_shape, activation="sigmoid"))

        adam = optimizers.Adam(learning_rate=1e-4,
                               clipnorm=0.5,
                               decay=(1e-4 / 100.0))
        model.compile(optimizer=adam,
                      loss="binary_crossentropy",
                      metrics=["accuracy"])
        self.model = model

    def train(self, X, y, sample_weights, n_epochs=100):
        callbacks = [
            keras.callbacks.EarlyStopping(monitor="val_loss", patience=3),
            keras.callbacks.History(),
        ]
        history = self.model.fit(
            x=X,
            y=y,
            epochs=n_epochs,
            shuffle=True,
            validation_split=0.2,
            batch_size=100,
            verbose=1,
            callbacks=callbacks,
            sample_weight=sample_weights,
        )
        return history

    def evaluate(self, X, y, batch_size, verbose):
        return self.model.evaluate(
            X, y, batch_size, verbose
        )

    def to_yaml(self):
        return self.model.to_yaml()
        
    def save(self, h5file):
        self.model.save(h5file)

    def predict(self, X, batch_size=100):
        return self.model.predict(X, batch_size)

    def error(self, X, y):
        return np.linalg.norm(y - self.model.predict(X))
