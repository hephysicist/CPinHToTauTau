from __future__ import annotations

import law
import math

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.util import maybe_import, InsertableDict
from columnflow.columnar_util import set_ak_column, flat_np_view, optional_column as optional

np = maybe_import("numpy")
ak = maybe_import("awkward")


logger = law.logger.get_logger(__name__)

@selector(
    uses={
        "PuppiMET.{covXX,covXY,covYY}",
    },
)
def met_cov_check(
    self: Selector,
    events: ak.Array,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    
    met = events.PuppiMET

    met_cov_mask = (
        np.isfinite(met.covXX) &
        np.isfinite(met.covXY) &
        np.isfinite(met.covYY) &
        ~np.isnan(met.covXX) &
        ~np.isnan(met.covXY) &
        ~np.isnan(met.covYY)
    )


    results = SelectionResult(
        steps={"met_cov_check": met_cov_mask},
    )

    return events, results