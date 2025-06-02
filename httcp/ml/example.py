# coding: utf-8

"""
Test model definition.
"""

from __future__ import annotations

import law
import order as od

from columnflow.types import Any
from columnflow.ml import MLModel
from columnflow.util import maybe_import, dev_sandbox
from columnflow.columnar_util import Route, set_ak_column

ak = maybe_import("awkward")
tf = maybe_import("tensorflow")

law.contrib.load("tensorflow")


class ExampleModel(MLModel):

    # mark the model as accepting only a single config
    single_config = True

    def setup(self):
        # dynamically add variables for the quantities produced by this model
        if f"{self.cls_name}.output" not in self.config_inst.variables:
            self.config_inst.add_variable(
                name=f"{self.cls_name}.output",
                null_value=-1,
                binning=(20, -1.0, 1.0),
                x_title=f"{self.cls_name} DNN output",
            )

    def sandbox(self, task: law.Task) -> str:
        return dev_sandbox("bash::$HTTCP_BASE/sandboxes/example.sh")

    def datasets(self, config_inst: od.Config) -> set[od.Dataset]:
        return {
            config_inst.get_dataset("st_tchannel_t_powheg"),
            config_inst.get_dataset("tt_sl_powheg"),
        }

    def uses(self, config_inst: od.Config) -> set[Route | str]:
        return {
            "Jet.pt", "Muon.pt",
        }

    def produces(self, config_inst: od.Config) -> set[Route | str]:
        return {
            f"{self.cls_name}.ouptut",
        }

    def output(self, task: law.Task) -> law.FileSystemDirectoryTarget:
        return task.target(f"mlmodel_f{task.branch}of{self.folds}", dir=True)

    def open_model(self, target: law.FileSystemDirectoryTarget) -> tf.keras.models.Model:
        return target.load(formatter="tf_keras_model")

    def train(
        self,
        task: law.Task,
        input: dict[str, list[dict[str, law.FileSystemFileTarget]]],
        output: law.FileSystemDirectoryTarget,
    ) -> None:
        # define a dummy NN
        x = tf.keras.Input(shape=(2,))
        a1 = tf.keras.layers.Dense(10, activation="elu")(x)
        y = tf.keras.layers.Dense(2, activation="softmax")(a1)
        model = tf.keras.Model(inputs=x, outputs=y)

        # the output is just a single directory target
        output.dump(model, formatter="tf_keras_model")

    def evaluate(
        self,
        task: law.Task,
        events: ak.Array,
        models: list[Any],
        fold_indices: ak.Array,
        events_used_in_training: bool = False,
    ) -> ak.Array:
        # fake evaluation
        events = set_ak_column(events, f"{self.cls_name}.output", 0.5)

        return events


# usable derivations
example = ExampleModel.derive("example", cls_dict={"folds": 2})



class HcpClassifierModel(MLModel):

    # mark the model as accepting only a single config
    single_config = True

    def setup(self):
        # dynamically add variables for the quantities produced by this model
        if f"{self.cls_name}.output" not in self.config_inst.variables:
            self.config_inst.add_variable(
                name=f"{self.cls_name}.output",
                null_value=-1,
                binning=(20, -1.0, 1.0),
                x_title=f"{self.cls_name} DNN output",
            )

    def sandbox(self, task: law.Task) -> str:
        return dev_sandbox("bash::$HTTCP_BASE/sandboxes/example.sh")

    def datasets(self, config_inst: od.Config) -> set[od.Dataset]:
        return {
            config_inst.get_dataset("h_ggf_htt_filtered"), #CM higgs signal
            config_inst.get_dataset("dy_lep_madgraph"),
        }

    def uses(self, config_inst: od.Config) -> set[Route | str]:
        ch_str = self.config_inst.channels.names()[0]
        return {
            f"hcand_{ch_str}.lep0.{var}"
            for var in ['pt','eta','phi']            
        } | {
            f"hcand_{ch_str}.lep1.{var}"
            for var in ['pt','eta','phi']            
        } | {
            f"hcand_{ch_str}.{var}"
            for var in ['delta_r','mt','pt','mass']            
        } | {
            f"PuppiMET.{var}"
            for var in ['pt','phi']
        }

    def produces(self, config_inst: od.Config) -> set[Route | str]:
        return {
            f"{self.cls_name}.ouptut",
        }

    def output(self, task: law.Task) -> law.FileSystemDirectoryTarget:
        return task.target(f"mlmodel_f{task.branch}of{self.folds}", dir=True)

    def build_model(self):
        # define a dummy NN
        x = tf.keras.Input(shape=(12,))
        a1 = tf.keras.layers.Dense(24, activation="elu")(x)
        a2 = tf.keras.layers.Dense(24, activation="elu")(a1)
        a3 = tf.keras.layers.Dense(12, activation="elu")(a2)
        y = tf.keras.layers.Dense(2, activation="softmax")(a3)
        model = tf.keras.Model(inputs=x, outputs=y)
        return model

    def open_model(self, target: law.FileSystemDirectoryTarget) -> tf.keras.models.Model:
        return target.load(formatter="tf_keras_model")

    def open_input_files(self, inputs):
        # contains files from all datasets
        events_of_datasets = inputs["events"][self.config_inst.name]

        # get datasets names
        # datasets = [dataset.label for dataset in self.datasets(self.config_inst)]

        # extract all columns from parquet files for all datasets and stack them
        all_events = []
        for dataset, parquet_file_targets in events_of_datasets.items():
            for parquet_file_target in parquet_file_targets:
                parquet_file_path = parquet_file_target["mlevents"].path
                events = ak.from_parquet(parquet_file_path)
                all_events.append(events)

        all_events = ak.concatenate(all_events)
        return all_events
    
    def prepare_events(self, events):
        # helper function to extract events and prepare them for training

        column_names = set(events.fields)
        input_features = set(self.input_features)
        target_features = set(self.target_features)

        # remove columns not used in training
        to_remove_columns = list(column_names - (input_features | target_features))

        for to_remove_column in to_remove_columns:
            print(f"removing column {to_remove_column}")
            events = remove_ak_column(events, to_remove_column)

        # ml model can't work with awkward arrays
        # we need to convert them to tf.tensors
        # this is done by following step chain:
        # ak.array -> change type to uniform type -> np.recarray -> np.array -> tf.tensor

        # change dtype to uniform type
        events = ak.values_astype(events, "float32")
        # split data in inputs and target
        input_columns = [events[input_column] for input_column in self.input_features]
        target_columns = [events[target_column] for target_column in self.target_features]

        # convert ak.array -> np.array -> bring in correct shape
        input_data = ak.concatenate(input_columns).to_numpy().reshape(
            len(self.input_features), -1).transpose()
        target_data = ak.concatenate(target_columns).to_numpy().reshape(
            len(self.target_features), -1).transpose()
        return tf.convert_to_tensor(input_data), tf.convert_to_tensor(target_data)
    
    
    def train(
        self,
        task: law.Task,
        input: dict[str, list[dict[str, law.FileSystemFileTarget]]],
        output: law.FileSystemDirectoryTarget,
    ) -> None:
        # define a dummy NN
        model = self.build_model()
        
        # get data tensors
        events = self.open_input_files(input)
        input_tensor, target_tensor = self.prepare_events(events)
        
        # setup everything needed for training
        optimizer = tf.keras.optimizers.SGD()
        model.compile(
            optimizer,
            loss="mse",
            steps_per_execution=10,
        )
        
        # train, throw model_history away
        _ = model.fit(
            input_tensor,
            target_tensor,
            epochs=5,
            steps_per_epoch=10,
            validation_split=0.25,
        )
        
        
        # the output is just a single directory target
        output.dump(model, formatter="tf_keras_model")

    def evaluate(
        self,
        task: law.Task,
        events: ak.Array,
        models: list[Any],
        fold_indices: ak.Array,
        events_used_in_training: bool = False,
    ) -> ak.Array:
        # fake evaluation
        events = set_ak_column(events, f"{self.cls_name}.output", 0.5)

        return events


# usable derivations
example = HcpClassifierModel.derive("example", cls_dict={"folds": 2})
