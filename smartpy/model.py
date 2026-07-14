import os
import json
import numpy as np
import pandas as pd

from ._initialise import initialise
from ._run import run
from ._finalise import finalise


def _assign_df_index(
        df: pd.DataFrame, var: str, regex: str, dtypes: list
) -> pd.DataFrame:
    # try to assign index if not already
    if isinstance(df.index, pd.RangeIndex):
        # try to guess names column
        if df.columns.dtype == 'str':
            # with headers, guess by label
            df_names = df.filter(regex=regex)
        elif df.columns.dtype == int:
            # without headers, guess by dtype
            df_names = df.select_dtypes(include=dtypes)
        else:
            raise RuntimeError(
                f"cannot guess {var!r} column, "
                f"please set {var!r} as index explicitly"
            )

        if df_names.size == 0:
            raise RuntimeError(
                f"{var!r} index not set "
                f"and cannot find {var} column"
            )
        elif df_names.columns.size > 1:
            raise RuntimeError(
                f"{var!r} index not set "
                f"and more than one {var!r} column found, "
                f"please set {var!r} as index explicitly"
            )

        df = df.set_index(df_names.columns[0])

    return df


def _resample_observations(df, idx) -> pd.DataFrame:
    # TODO: functional but a bit slow to run, find a faster approach

    # target intervals (assumed regular)
    idx_interval = idx[1] - idx[0]
    idx_starts = idx - idx_interval

    # current intervals (may be irregular)
    df_ends = df.index.values
    df_starts = df.index.to_series().shift()
    # /!\ assume first interval is same as second interval
    df_starts.iloc[0] = df.index[0] - (df.index[1] - df.index[0])
    df_starts = df_starts.values

    arr = df.values

    def weighted_avg(s, e):
        # only keep df intervals that can overlap
        i = np.searchsorted(df_ends, s, side="right")
        j = np.searchsorted(df_starts, e, side="left")

        # compute overlaps in seconds
        overlaps = (
            # overlap starts at the later start time
            np.minimum(df_ends[i:j], e)
            # overlap ends at the earlier end time
            - np.maximum(df_starts[i:j], s)
        ).astype("timedelta64[s]").astype(float)

        # use overlaps as weights to compute average
        return np.average(arr[i:j], axis=0, weights=overlaps)

    return pd.DataFrame(
        [weighted_avg(s, e) for s, e in zip(idx_starts, idx)],
        index=idx,
        columns=df.columns,
    )


def _process_df_timeseries(
        df: pd.DataFrame, var: str, target_index: pd.DatetimeIndex = None
) -> pd.DataFrame:
    # check type
    if not isinstance(df, pd.DataFrame):
        raise TypeError(
            f"{var!r} data must be provided as pandas.DataFrame"
        )

    # assign temporal info as index if not already
    df = _assign_df_index(
        df, var, '[Dd][Aa][Tt][Ee]|[Tt][Ii][Mm][Ee]', [str, np.datetime64]
    )

    # make sure it is a datetime index
    df.index = pd.to_datetime(df.index)
    df = df.sort_index()

    # in case of observed data, resample onto forcing index
    if target_index is not None:
        df = _resample_observations(df, target_index)

    # convert [kg m-2 timedelta-1] into [kg m-2 s-1]
    deltas = df.index.diff()[1:]
    if deltas.min() != deltas.max():
        raise RuntimeError(
            f"{var!r} data must feature a constant temporal frequency"
        )

    if target_index is None:
        df = df.asfreq(deltas[0])
        df = df / deltas[0].total_seconds()

    return df


def _process_df_parameters(df: pd.DataFrame) -> pd.DataFrame:
    # check type
    if not isinstance(df, pd.DataFrame):
        raise TypeError(
            "parameters must be provided as pandas.DataFrame"
        )

    # assign parameter names as index if not already
    df = _assign_df_index(df, 'parameters', '[Nn][Aa][Mm][Ee]', [str])

    # make sure it contains all parameters
    names = np.array(['T', 'C', 'H', 'D', 'S', 'Z', 'SK', 'FK', 'GK', 'RK'])
    found = df.index.isin(names)

    if not found.sum() == len(names):
        raise KeyError(
            f"{names[~found]!r} parameters not found in index"
        )

    # make sure order is correct
    df = df.reindex(names)

    return df


def _get_number_basins(**kwargs) -> int:
    # get number of columns
    n_cols = {var: df.shape[1] for var, df in kwargs.items()}

    if len(set(n_cols.values())) == 1:
        return list(n_cols.values())[0]
    else:
        raise RuntimeError(
            f"number of basins not matching: {n_cols!r}"
        )


class Model(object):

    app_path = os.path.dirname(__file__)

    with open(os.path.join(app_path, 'model.json'), 'r') as f:
        _meta = json.load(f)

    def __init__(
            self,
            rain: pd.DataFrame, pet: pd.DataFrame, area: np.ndarray | float,
            flow: pd.DataFrame = None
    ):
        # assign forcing data time series
        self.rain = rain
        self.pet = pet
        if not self.rain.index.equals(self.pet.index):
            raise RuntimeError("'rain' and 'pet' time index are not equal")
        else:
            self.timedelta = (
                    self.rain.index[1] - self.rain.index[0]
            ).total_seconds()

        # assign basin(s) area(s)
        self.area = area

        # assign observed data time series
        self.flow = flow

        # determine number of basins
        self.nx = _get_number_basins(
            rain=self.rain, pet=self.pet, area=pd.DataFrame([self.area]),
            **(dict(flow=self.flow) if flow is not None else dict())
        )

    @property
    def rain(self):
        """Return the rainfall timeseries as a `pandas.DataFrame`."""
        return self._rain.copy()

    @rain.setter
    def rain(self, rain):
        self._rain = _process_df_timeseries(rain, 'rain')

    @property
    def pet(self):
        """Return the potential evapotranspiration timeseries as a `pandas.DataFrame`."""
        return self._pet.copy()

    @pet.setter
    def pet(self, pet):
        self._pet = _process_df_timeseries(pet, 'pet')

    @property
    def area(self):
        """Return the basin(s) area(s) as a `numpy.array`."""
        return self._area.copy()

    @area.setter
    def area(self, area):
        self._area = np.asarray(area)

    @property
    def flow(self):
        """Return the streamflow timeseries as a `pandas.DataFrame`."""
        return self._flow.copy()

    @flow.setter
    def flow(self, flow):
        if flow is not None:
            self._flow = _process_df_timeseries(
                flow, 'flow', target_index=self.rain.index
            )
        else:
            self._flow = None

    def simulate(
            self, parameters: pd.DataFrame,
            start: str | pd.Timestamp = None, end: str | pd.Timestamp = None,
            istart: str | pd.Timestamp = None, iend: str | pd.Timestamp = None,
            _spinup: bool = False
    ):
        dtype = np.float64

        # gather inputs for relevant period
        inputs = {
            'rainfall_flux':
                self.rain.loc[start:end, :].values,
            'potential_evapotranspiration_flux':
                self.pet.loc[start:end, :].values
        }
        nt = inputs['rainfall_flux'].shape[0]

        # allocate memory for states
        states = {
            name: np.zeros(
                (nt + 1, self.nx) if 'divisions' not in attrs else
                (nt + 1, self.nx, attrs['divisions']),
                dtype=dtype
            ) for name, attrs in self._meta['states'].items()
        }

        # spin up to set initial conditions
        if istart is not None and iend is not None:
            istates = self.simulate(
                parameters, istart, iend, _spinup=True
            )
            for name, attrs in self._meta['states'].items():
                states[name][0, ...] = istates[name][-1, ...]

        # gather parameter values
        parameters = _process_df_parameters(parameters)
        _ = _get_number_basins(rain=self.rain, parameters=parameters)
        parameters = {
            name: parameters.loc[attrs['name']].values
            for name, attrs in self._meta['parameters'].items()
        }

        # gather constants values
        constants = {
            'timedelta':
                self.timedelta,
            'drainage_area':
                self.area,
            'rho_water':
                1e3
        }

        # allocate memory for outputs
        outputs = {
            name: np.zeros((nt, self.nx), dtype=dtype)
            for name, attrs in self._meta['outputs'].items()
        }

        # create placeholders for internals
        internals = {
            name: np.zeros((1, self.nx), dtype=dtype)
            for name, attrs in self._meta['internals'].items()
        }

        # call initialise/run/finalise functions
        initialise(
            **parameters, **states, **constants
        )
        run(
            **inputs, **parameters, **states, **constants,
            **internals, **outputs, nt=nt
        )
        finalise()

        return {
            name: outputs[name] for name, attrs in self._meta['outputs'].items()
        } if not _spinup else {
            name: states[name] for name, attrs in self._meta['states'].items()
        }


if __name__ == '__main__':
    df_rain = pd.read_csv('../tests/data/in/Catchment/Catchment.rain')
    df_rain = _process_df_timeseries(df_rain, 'rain')

    df_pet = pd.read_csv('../tests/data/in/Catchment/Catchment.peva')
    df_pet = _process_df_timeseries(df_pet, 'pet')

    df_flow = pd.read_csv('../tests/data/in/Catchment/Catchment.flow', index_col=0)
    # df_flow = _process_df_timeseries(df_flow, 'flow', target_index=df_rain.index)

    # # single basin
    # m = Model(
    #     rain=df_rain,
    #     pet=df_pet,
    #     area=1,
    #     flow=df_flow
    # )

    # multiple basins
    m = Model(
        rain=pd.concat([df_rain] * 5, ignore_index=True, axis=1),
        pet=pd.concat([df_pet] * 5, ignore_index=True, axis=1),
        area=[100] * 5,
        # flow=pd.concat([df_flow] * 5, ignore_index=True, axis=1)
    )

    # TODO: check that a Monte Carlo simulation works with broadcasting
    #       i.e. that by providing several parameters columns (sets)
    #       it can produce as many outputs columns without having to
    #       provide several inputs (forcings/observations) columns
    df_params = pd.read_csv('../tests/data/in/Catchment/Catchment.parameters', index_col=0)

    r = m.simulate(
        # my_df_params,
        pd.concat([df_params] * 5, ignore_index=True, axis=1),
        start='2003-01-01', end='2006-12-31',
        istart='2000-01-01', iend='2002-12-31'
    )

    print(r['river_discharge_flux'])
