# coding: utf-8

import os
import time
import pathlib
from multiprocessing import Process, Pipe
from multiprocessing.connection import Connection
from dataclasses import dataclass
from typing import Any

import law
import order as od

from columnflow.types import Any
from columnflow.ml import MLModel
from columnflow.util import maybe_import, dev_sandbox,InsertableDict
from columnflow.columnar_util import Route, set_ak_column
ak = maybe_import("awkward")
xgb = maybe_import("xgboost")

STOP_SIGNAL = "STOP"


class XGBEvaluator:
    """
    TensorFlow model evaluator that runs in a separate process with support for multiple models.

    .. code-block:: python

        evaluator = TFEvaluator()
        evaluator.add_model("model_name", "path/to/model")
        with evaluator:
            result = evaluator("model_name", input_data)
    """

    @dataclass
    class Model:
        name: str
        path: str
        pipe: Connection = None
        signature_key: str = ""

    def __init__(self) -> None:
        super().__init__()

        self._models: dict[str, XGBEvaluator.Model] = {}
        self._p: Process | None = None

        self.delay = 0.2
        self.silent = False

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
        self.stop()

    def __del__(self) -> None:
        self.stop()

    def __call__(self, *args, **kwargs) -> Any:
        return self.evaluate(*args, **kwargs)

    @property
    def running(self) -> bool:
        return self._p is not None

    def add_model(self, name: str, path: str, signature_key = "") -> None:
        if self.running:
            raise ValueError("cannot add models while running")
        if name in self._models:
            raise ValueError(f"model with name '{name}' already exists")

        # normalize path
        path = str(path)
        path = os.path.expandvars(os.path.expanduser(path))
        path = os.path.abspath(os.path.abspath(path))

        # add it
        self._models[name] = XGBEvaluator.Model(name=name, path=path, signature_key=signature_key)

    def start(self) -> None:
        if self.running:
            raise ValueError("process already started")

        # build the subprocess config
        config = []
        for model in self._models.values():
            parent_pipe, child_pipe = Pipe()
            model.pipe = parent_pipe
            config.append({"name": model.name, "path": model.path, "pipe": child_pipe})

        # create and start the process
        self._p = Process(
            target=_xgb_evaluate,
            args=(config,),
            kwargs={"delay": self.delay, "silent": self.silent},
        )
        self._p.start()

    def evaluate(self, name: str, *args, **kwargs) -> Any:
        if not self.running:
            raise ValueError("process not started")

        # get the model
        if name not in self._models:
            raise ValueError(f"model with name '{name}' does not exist")
        model = self._models[name]

        # evaluate and send back result
        model.pipe.send((args, kwargs))  # type: ignore[union-attr]
        return model.pipe.recv()  # type: ignore[union-attr]

    def stop(self, timeout = 5) -> None:
        # stop and remove model pipes
        for model in self._models.values():
            if model.pipe is not None:
                model.pipe.send(STOP_SIGNAL)
                model.pipe.close()
                model.pipe = None

        # nothing to do when not running
        if not self.running:
            return

        # join to wait for normal termination
        if self._p.is_alive():
            self._p.join(timeout)

            # kill if still alive
            if self._p.is_alive():
                self._p.kill()

        # reset
        self._p = None


def _xgb_evaluate(
    config: list[dict[str, Any]],
    /,
    *,
    delay = 0.2,
    silent: bool = False,
) -> None:
    _print = (lambda *args, **kwargs: None) if silent else print

    _print("importing xgboost ...")
    import numpy as np
    import xgboost as xgb  # type: ignore[import-not-found,import-untyped]
    _print("done")

    @dataclass
    class Model:
        name: str
        path: str
        pipe: Connection
        signature_key: str = ""
        model: Any = None

        @classmethod
        def new(cls, config: dict[str, Any]):
            for attr in ("name", "path", "pipe"):
                if attr not in config:
                    raise ValueError(f"missing field '{attr}' in model config")
            if not os.path.exists(config["path"]):
                raise FileNotFoundError(f"model file '{config['path']}' does not exist")
            if not isinstance(config["pipe"], Connection):
                raise TypeError(f"'pipe' {config['pipe']} not of type '{Connection}'")
            return cls(
                name=config["name"],
                path=config["path"],
                pipe=config["pipe"],
                signature_key=config.get("signature_key", ""),
            )

        def load(self) -> None:
            sig_msg = f" (signature '{self.signature_key}')" if self.signature_key else ""
            _print(f"loading model '{self.name}'{sig_msg} from {self.path} ...")
            self.model = xgb.XGBClassifier()
            self.model.load_model(self.path)
            _print("done")

        def evaluate(self, *args, **kwargs) -> np.ndarray:
            return self.model.predict_proba(*args, **kwargs)

        def clear(self) -> None:
            _print(f"clearing model '{self.name}'")
            self.model = None
            self.pipe.close()

    # convert to model objects
    models = [Model.new(item) for item in config]

    # load model objects
    for model in models:
        model.load()

    # helper for gracefully shutting down
    def shutdown() -> None:
        for model in models:
            model.clear()
        models.clear()

    # start loop listening for data
    while models:
        remove_models: list[int] = []
        for i, model in enumerate(models):
            # skip if there is no data to process
            if not model.pipe.poll():
                continue

            # get data and process
            data = model.pipe.recv()
            if isinstance(data, tuple) and len(data) == 2:
                # evaluate
                try:
                    args, kwargs = data
                    result = model.evaluate(*args, **kwargs)
                except:
                    shutdown()
                    raise
                # send back result
                model.pipe.send(result)

            elif data == STOP_SIGNAL:
                # remove model
                model.clear()
                remove_models.append(i)

            else:
                raise ValueError(f"unexpected data type {type(data)}")

        # reduce models and sleep
        models = [model for i, model in enumerate(models) if i not in remove_models]
        time.sleep(delay)
        





















# class XGBEvaluator(MLModel):

#     # mark the model as accepting only a single config
#     single_config = True
    
#     def __init__(self) -> None:
#         super().__init__()

#         self._models: dict[str, XGBEvaluator.Model] = {}
#         self._p: Process | None = None

#         self.delay = 0.2
#         self.silent = False
    
#     def __enter__(self) -> XGBEvaluator:
#         self.start()
#         return self

#     def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
#         self.stop()

#     def __del__(self) -> None:
#         self.stop()
        
#     def __call__(self, *args, **kwargs) -> Any:
#         return self.evaluate(*args, **kwargs)
        
#     @property
#     def running(self) -> bool:
#         return self._p is not None
    
#     def start(self) -> None:
#         if self.running:
#             raise ValueError("process already started")

#         # build the subprocess config
#         config = []
#         for model in self._models.values():
#             parent_pipe, child_pipe = Pipe()
#             model.pipe = parent_pipe
#             config.append({"name": model.name, "path": model.path, "pipe": child_pipe})

#         # create and start the process
#         self._p = Process(
#             target=_tf_evaluate,
#             args=(config,),
#             kwargs={"delay": self.delay, "silent": self.silent},
#         )
#         self._p.start()

#     def evaluate(self, name: str, *args, **kwargs) -> Any:
#         if not self.running:
#             raise ValueError("process not started")

#         # get the model
#         if name not in self._models:
#             raise ValueError(f"model with name '{name}' does not exist")
#         model = self._models[name]

#         # evaluate and send back result
#         model.pipe.send((args, kwargs))  # type: ignore[union-attr]
#         return model.pipe.recv()  # type: ignore[union-attr]

#     def stop(self, timeout: int | float = 5) -> None:
#         # stop and remove model pipes
#         for model in self._models.values():
#             if model.pipe is not None:
#                 model.pipe.send(STOP_SIGNAL)
#                 model.pipe.close()
#                 model.pipe = None

#         # nothing to do when not running
#         if not self.running:
#             return

#         # join to wait for normal termination
#         if self._p.is_alive():
#             self._p.join(timeout)

#             # kill if still alive
#             if self._p.is_alive():
#                 self._p.kill()

#         # reset
#         self._p = None
    
    
    

#     # def setup(self):
#     #     # dynamically add variables for the quantities produced by this model
#     #     if f"{self.cls_name}.output" not in self.config_inst.variables:
#     #         self.config_inst.add_variable(
#     #             name=f"{self.cls_name}.output",
#     #             null_value=-1,
#     #             binning=(20, -1.0, 1.0),
#     #             x_title=f"{self.cls_name} DNN output",
#     #         )

#     # def sandbox(self, task: law.Task) -> str:
#     #     return dev_sandbox("bash::$HTTCP_BASE/sandboxes/example.sh")

#     # def datasets(self, config_inst: od.Config) -> set[od.Dataset]:
#     #     return {
#     #         config_inst.get_dataset("st_tchannel_t_powheg"),
#     #         config_inst.get_dataset("tt_sl_powheg"),
#     #     }

#     # def uses(self, config_inst: od.Config) -> set[Route | str]:
#     #     return {
#     #         "Jet.pt", "Muon.pt",
#     #     }

#     # def produces(self, config_inst: od.Config) -> set[Route | str]:
#     #     return {
#     #         f"{self.cls_name}.ouptut",
#     #     }

#     # def output(self, task: law.Task) -> law.FileSystemDirectoryTarget:
#     #     return task.target(f"mlmodel_f{task.branch}of{self.folds}", dir=True)

#     # def open_model(self, target: law.FileSystemDirectoryTarget) -> tf.keras.models.Model:
#     #     return target.load(formatter="tf_keras_model")

#     # def train(
#     #     self,
#     #     task: law.Task,
#     #     input: dict[str, list[dict[str, law.FileSystemFileTarget]]],
#     #     output: law.FileSystemDirectoryTarget,
#     # ) -> None:
#     #     # define a dummy NN
#     #     x = tf.keras.Input(shape=(2,))
#     #     a1 = tf.keras.layers.Dense(10, activation="elu")(x)
#     #     y = tf.keras.layers.Dense(2, activation="softmax")(a1)
#     #     model = tf.keras.Model(inputs=x, outputs=y)

#     #     # the output is just a single directory target
#     #     output.dump(model, formatter="tf_keras_model")

#     # def evaluate(
#     #     self,
#     #     task: law.Task,
#     #     events: ak.Array,
#     #     models: list[Any],
#     #     fold_indices: ak.Array,
#     #     events_used_in_training: bool = False,
#     # ) -> ak.Array:
#     #     # fake evaluation
#     #     events = set_ak_column(events, f"{self.cls_name}.output", 0.5)

#     #     return events


# # usable derivations
# #example = ExampleModel.derive("example", cls_dict={"folds": 2})
