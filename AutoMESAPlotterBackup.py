# Author:  Jess Vriesema, Nathan Steenwyk and Ethan Webber
# Date:    4 June 2026
# Purpose: Define plotting tools for use with AutoMESA.


import numpy as np
import matplotlib.pyplot as plt
import mesa_reader as mr
import copy
from itertools import product as iter_product
from MESA_models import MESAModelGrid
from MESA_models import MESAModel
from dataclasses import dataclass
from typing import Callable
from numpy.typing import NDArray
from sympy import *



"""
------------------------------------------------------------------------
DESIGN CONSIDERATIONS
------------------------------------------------------------------------

Class Design/Structure:
----------------------
- Encapsulate data as much as possible, so that the user does not have to 
    worry about how the data is stored or accessed. Look up tips for 
    implementing encapsulation in Python. 

Getting Data:
------------
- We want to be able to easily get the data from the MESA output files for
    any model in the model grid, without having to worry about how the data 
    is stored or accessed. We want to be able to easily get the data for any
    variable from any MESA output file for any model in the model grid.

    

    def get_profile_data_with_auto_mesa(self, 
                 gridVars: tuple(tuple[str,float]),
                 profile:str, 
                 variables: str|list[str]=None, 
                 zones: int|list[int]=None|None) -> np.ndarray:


Examples of loading data from MESA output files:
-----------------------------------------------
    # Grabs temp and pressure data from profile1.data for the model with 
    #   mdot=1.0 and Ltrans=0.02 in the model grid:
    data = plotter.get_profile_data_with_auto_mesa( 
                    gridVars=(("mdot", 1.0), ("Ltrans", 0.02)),
                    profile="profile_600.data",
                    variables=["T", "P"],
                    zones=None )
    # Plot P vs T for that model:
    plt.plot( data["T"], data["P"] )

    # Grabs the surface temp and pressure for multiple models:
    data = plotter.get_profile_data_with_auto_mesa(
                    gridVars=(("mdot", None), ("Ltrans", 0.02)),
                    profile="profile_600.data",
                    variables=["T", "P"],
                    zones=None )
    # Plot Psurf vs Tsurf for each model on the same plot for different mdot values:
    for mdot in data.keys():
        plt.plot( data[mdot]["T"], data[mdot]["P"] )

    # Grabs the star's outer radius at the last profile for a whole 2D grid of models:
    data = plotter.get_profile_data_with_auto_mesa(
                    gridVars=(("mdot", None), ("Ltrans", None)),
                    profile="profile_end.data",
                    variables=["R"],
                    zones=1 )
    # Plot Rsurf vs mdot with lines for each Ltrans value:
    # TODO: this syntax needs work
    for l_trans in data.keys():
        plt.plot( data[("Ltrans",l_trans)].keys("mdot"), 
                  data[(("Ltrans",l_trans),("mdot",None))],
                  'o-', label=f"Ltrans = {l_trans}")

        

Data Handling:
-------------
- Lazy evaluation of data????  -- NOT NOW (or not yet?)
    - Only load the data from MESA output files when it is actually needed 
        for plotting. Do NOT load all the data from all the MESA output 
        files at once, since this could be a huge amount of data and 
        could cause memory issues later on with bigger data sets. 
    - Since there will potentially be many MESA models each with potentially 
        many output files, we want to be careful about how much data we load 
        into memory at once.
    - Bonus points if there could be a setting that enables lazy evaluation 
        of data, so that the user can choose whether to load all the data at 
        once or not. Why? Because lazy evaluation can be slower than loading 
        all the data at once, so for smaller data sets it might be better to 
        just load all the data at once.

Plotting:
--------
- We want to be able to easily plot any quantity from the MESA output files
    against any other quantity from the MESA output files.
- We want to be able to easily plot multiple models on the same plot, for
    comparison.
- We want to be able to easily customize the plots (e.g. labels, colors, etc.).
- We want to be able to do (scatter plots, line plots) as well as 
    (pcolor plots, contour plots, etc.). Both have different inputs. 
"""



# Note: frozen=True makes the dataclass immutable, which is good for our use 
#      case since we want the records to be read-only after they are created.
@dataclass(frozen=True)
class MultiDimSliceRecord:
    """One row from a slice query: full coordinate key and stored value."""

    key: tuple[tuple[str, object], ...]
    value: object


class MultiDimSliceResult:
    """List-like slice output with helper access to varying key values."""

    __slots__ = ["_records"]

    def __init__(self, records: list[MultiDimSliceRecord]) -> None:
        self._records = records

    def __iter__(self):
        for rec in self._records:
            yield rec.value

    def __len__(self) -> int:
        return len(self._records)

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            return MultiDimSliceResult(self._records[idx])
        return self._records[idx].value

    def values(self) -> list[object]:
        """Return the sliced data values as a plain list."""
        return [rec.value for rec in self._records]

    def records(self) -> list[MultiDimSliceRecord]:
        """Return full records including coordinates and values."""
        return list(self._records)

    def axis_values(self, key_name: str) -> list[object]:
        """Return all values for one key (e.g., key2) across the slice."""
        result = []
        for rec in self._records:
            rec_map = dict(rec.key)
            if key_name in rec_map:
                result.append(rec_map[key_name])
        return result


class MultiDimDict:
    """Dictionary-like container keyed by N dimensions.

    Keys are tuples of `(name, value)` pairs, for example:
        data[("key1", 1), ("key2", "A")] = payload

    Slicing is supported by passing `:` (`slice(None)`) in one or more
    dimensions, for example:
        out = data[("key1", 1), :]
        out.values()         # sliced data list
        out.axis_values("key2")  # corresponding varying key2 values
    """

    __slots__ = ["_store", "_key_order"]

    def __init__(self, key_order: list[str] | tuple[str, ...] | None = None) -> None:
        self._store: dict[tuple[tuple[str, object], ...], object] = {}
        self._key_order = list(key_order) if key_order is not None else None

    def __len__(self) -> int:
        return len(self._store)

    def __contains__(self, key) -> bool:
        norm = self._normalize_query(key)
        return norm in self._store

    def __setitem__(self, key, value) -> None:
        norm = self._normalize_exact_key(key)
        self._store[norm] = value

    def __getitem__(self, query):
        norm = self._normalize_query(query)
        if self._is_exact_query(norm):
            return self._store[norm]

        records: list[MultiDimSliceRecord] = []
        for stored_key, stored_val in self._store.items():
            if self._matches(norm, stored_key):
                records.append(MultiDimSliceRecord(stored_key, stored_val))

        if not records:
            raise KeyError(f"No entries matched query: {query}")

        return MultiDimSliceResult(records)

    def get(self, query, default=None):
        try:
            return self[query]
        except KeyError:
            return default

    def keys(self):
        return self._store.keys()

    def values(self):
        return self._store.values()

    def items(self):
        return self._store.items()

    def clear(self) -> None:
        self._store.clear()

    def key_values(self, key_name: str) -> list[object]:
        """Return unique values seen for one key dimension."""
        vals = []
        for key in self._store.keys():
            mapping = dict(key)
            if key_name in mapping and mapping[key_name] not in vals:
                vals.append(mapping[key_name])
        return vals

    def _normalize_exact_key(self, key) -> tuple[tuple[str, object], ...]:
        parsed = self._parse_query(key)
        for _, val in parsed:
            if isinstance(val, slice):
                raise TypeError("Slice is not allowed when setting a value.")

        return self._ordered_key(parsed)

    def _normalize_query(self, query) -> tuple[tuple[str, object], ...]:
        parsed = self._parse_query(query)
        return self._ordered_key(parsed)

    def _parse_query(self, query) -> list[tuple[str, object]]:
        if not isinstance(query, tuple):
            query = (query,)

        parsed = []
        for item in query:
            if item == slice(None):
                parsed.append(("*", slice(None)))
                continue

            if isinstance(item, tuple) and len(item) == 2 and isinstance(item[0], str):
                parsed.append((item[0], item[1]))
                continue

            raise TypeError(
                "Each index must be ':' or a (key_name, value) pair; "
                f"got {item!r}."
            )

        return parsed

    def _ordered_key(self, parsed: list[tuple[str, object]]) -> tuple[tuple[str, object], ...]:
        wildcard_count = sum(1 for name, _ in parsed if name == "*")

        names = [name for name, _ in parsed if name != "*"]
        if len(names) != len(set(names)):
            raise ValueError("Duplicate key names in index query are not allowed.")

        if self._key_order is None:
            if wildcard_count > 0:
                raise ValueError("Cannot use ':' before key order is known.")
            self._key_order = [name for name, _ in parsed]

        if len(parsed) != len(self._key_order):
            raise ValueError(
                f"Expected {len(self._key_order)} dimensions, got {len(parsed)}."
            )

        by_name = {name: val for name, val in parsed if name != "*"}
        wildcard_positions = [i for i, (name, _) in enumerate(parsed) if name == "*"]

        if wildcard_positions:
            missing = [name for name in self._key_order if name not in by_name]
            if len(missing) != len(wildcard_positions):
                raise ValueError(
                    "Could not map ':' positions to key names; provide all key names "
                    "except those intended as wildcards."
                )
            for name in missing:
                by_name[name] = slice(None)

        if set(by_name.keys()) != set(self._key_order):
            raise ValueError(
                f"Index names must match key_order={self._key_order}; got {sorted(by_name.keys())}."
            )

        return tuple((name, by_name[name]) for name in self._key_order)

    @staticmethod
    def _is_exact_query(query: tuple[tuple[str, object], ...]) -> bool:
        return all(not isinstance(val, slice) for _, val in query)

    @staticmethod
    def _matches(query, stored_key) -> bool:
        for (_, qval), (_, sval) in zip(query, stored_key):
            if isinstance(qval, slice):
                continue
            if qval != sval:
                return False
        return True










class AutoMESAPlotter:
    # Instance variables:
    __slots__ = ["modelGrid", "runName", "localAutoMESADir"]  # TODO: Add more instance variables as needed
    

    # For Nathan
    def __init__( self, localAutoMESADir:str, runName:str, AutoMESA_gridfile:str) -> None:
        # TODO: Fill this in
        self.modelGrid = MESAModelGrid(AutoMESA_gridfile, runName, localAutoMESADir)
        self.runName = runName
        self.localAutoMESADir = localAutoMESADir
        return
    

    def load_data( self, datafile:str ) -> None:
        """Only loads the data from a MESA output file that is actually needed."""
        # TODO: Fill this in
        return
    

    # For Nathan
    def get_model(self, gridVarDict) -> MESAModel: 
        """Returns the MESAModel object associated with this AutoMESAPlotter.
        Suppose you have a MESAModelGrid with grid variables "initial_mass" 
        and "initial_metallicity". Then you could call this method with 
        gridVarDict={"initial_mass": 1.0, "initial_metallicity": 0.02} to get 
        the MESAModel object associated with the model that has initial mass 
        of 1.0 and initial metallicity of 0.02. If no model in the model 
        grid has those grid variable values, then this method should raise 
        an error."""
        
            
        modelIndex = []
        for key, value in gridVarDict.items():
            if key not in self.modelGrid.varNames:
                found = False
                for nicknamed_key, nickname in self.modelGrid.nicknamedVars.items():
                    if key == nickname:
                        key = nicknamed_key
                        found = True
                        break
                if not found:
                    raise TypeError(f"gridVars name {key} is invalid")
            
            i = self.modelGrid.varNames.index(key)
            arr = self.modelGrid.varValues[i]
            index = next((i for i, val in enumerate(arr) if val == value), None)
            
            if index is None:
                raise ValueError(f"Value {value} not found under variable name {key}.")
            
            modelIndex.append(index)

        return self.modelGrid.models[tuple(modelIndex)]




    # For Ethan
    def get_data_simple(self, 
                 autoMESADir:str, 
                 modelName:str, 
                 profile:str, 
                 variables: str|list[str]=None, 
                 zones: int|list[int]|None=None) -> np.ndarray:
        """Returns the data series at all levels for all of the requested 
        variables from a given MESA output file.

        PARAMETERS:
        
        autoMESADir (str): 
        The AutoMESA directory in which the MESA models are located.
        
        modelName (str):
        The name of the MESA model. This is needed to find the MESA output file.

        variables (str or list of str): 
            The variable(s) to get the data for. 
            If variables is a string, then return the data series for that variable.
            If variables is a list of strings, then return a dictionary mapping 
              each variable to its data.
        
        zones (int, list of ints, or None): 
           The zones of the data series to return. If zones is None, 
           then return everything. If it is just an integer, return just 
           one value. If zones is a list of integers, then return only the 
           data at those zones. All zones should be non-negative integers. 
           If any index is out of bounds for the data series, then raise 
           an error.
        """

        # 1. Check that input parameters are of the correct type.
        #    (We could also check that the model and profile actually exist, 
        #    but that might be too much for now.)
        if not isinstance(autoMESADir, str):
            raise TypeError(f"autoMESADir must be a string, but was a {type(autoMESADir)} instead.")
        if not isinstance(modelName, str):
            raise TypeError(f"modelName must be a string, but was a {type(modelName)} instead.")
        if not isinstance(profile, str):
            raise TypeError(f"profile must be a string, but was a {type(profile)} instead.")
        if not isinstance(zones, (int, list, type(None))):
            raise TypeError(f"zones must be an int, a list of ints, or None, but was a {type(zones)} instead.")

        # 2. Use mesa_reader to read in the data from the specified MESA  
        #    output file.
        md = mr.MesaData(f'{autoMESADir}/{modelName}/LOGS/{profile}')
        
        
        if isinstance(variables, list):
            #Create an empty dicionary to hold the extracted arrays for each variable
            extracted_data = {}
            var_list = [variables] if isinstance(variables, str) else variables
            for var in var_list:
                #The in_data method checks if the variable is in the MESA output file
                if not md.in_data(var): 
                    raise ValueError(f"Variable '{var}' not found in MESA output file '{profile}' for model '{modelName}' in AutoMESA directory '{autoMESADir}'.")
                else:
                    data_array = md.data(var)
                    #If zones is not None, then return only the data at those specified zones.
                    if isinstance(zones, list) or zones is None:
                        # We need to return a list of values at those zones
                        extracted_data[var] = data_array[zones] if zones is not None else data_array
                    elif isinstance(zones, int):
                        # We need to return just one value at that index
                         extracted_data[var] = data_array[zones] 
                    else:
                        # raise an error! zones should be an int, a list of ints, or None.
                        raise TypeError(f"zones must be an int, a list of ints, or None, but was a {type(zones)} instead.")
    
            #Returns a dictionary that maps each variable to its data. 
            return extracted_data 
        elif isinstance(variables, str):
             #Checks if the single string exsits in the MESA file output file. 
            if not md.in_data(variables):
                raise ValueError(f"Variable '{variables}' not found in MESA output file '{profile}' for model '{modelName}' in AutoMESA directory '{autoMESADir}'.")
            else:
                data_array = md.data(variables)
                #If zones is not None, then return only the data at those specified zones.
                if zones is not None: 
                    data_array = data_array[zones]
                return data_array
        else:
            raise TypeError(f"variables must be a string or a list of strings, but was a {type(variables)} instead.")
        
            
        # 3. Return the data series for the specified variable.
        #    If variable is just a single string, then return a np.ndarray
        #       with the data for that variable.
        
        #    If the variable(s) requested do not exist in the MESA output file, then
        #       raise an error.

        ...  # placeholder that tells us we aren't done with the method yet!
    


    # For Nathan
    def get_profile_data_with_auto_mesa(self, 
                 gridVars: dict[str, float|list[float]|None],
                 profile: str, 
                 variables: str|list[str]|None=None, 
                 zones: int|list[int]|list[str]|Callable[[mr.MesaData], NDArray[np.bool_]]|None=None,
                 simplify: bool=True) -> float|int|np.ndarray[float]|dict:
        """
        Returns the data series at specified levels for given variables from 
        given MESA models.

        
        Parameters
        ----------
        gridVars : dict[str, float|list[float]|None]
            Dictionary defining a desired sub-grid of models by mapping variables to values.
            Any nicknames defined in ``self.modelGrid.nicknamedVars` may be used in place of 
            the original names. Tied variables are also accepted.
        
        profile : str
            A string with the name of the desired profile file.

        variables : str, list[str], or None
            The variable(s) to retrieve. If variables is None, then all variables
            in the profile are returned.
        
        zones : int, list[int], list[str], Callable[[mr.MesaData], NDArray[np.bool_]], or None, default=None
            The zones of the data series to return. If zones is None, then all zones are 
            returned. Zones can take a list of strings giving an equation and its variables
            to determine which zones are output. Can also take a callable to determine which
            zones are output via a boolean mask.
        
        simplify : bool, default=True
            If true, data structure returned is simplified to remove unneccessary nesting.
            See Notes for more.
        
            
        Returns
        -------
        float, int, np.ndarray[float], or dict
            The data of the specified models, variables, and zones from the specified profile file.
            The exact return type and its structure depends on whether one or multiple of each 
            parameter is specified.
            
        
        Notes
        -----
        The full return structure is arranged hierarchically as follows:
            * dictionary with data from all models, keys specify model
                * dictionary with data from one model, keys specify variable
                    * array with values of one variable from model

        The data of individual models is keyed by nested tuples of (gridVarName, value) pairs.
        The order the variables is the order of insertion of the same variables into the
        dictionary given to gridVars. If a nickname or tied variable name is given, those will
        be used instead of the grid variable.

        When simplify=True, any level containing only a single item is removed. This corresponds
        to whether one or multiple values are specified for the parameters. When simplify=False,
        the full data structure is always returned.

        To give an equation to parameter zones, it must be in the form [equation, var1, var2, etc.].
        The variable names given in the list can be any string, and should be used in the equation. 
        The equation is parsed by sympy.parse_expr(), so the syntax is mostly standard python syntax.
        Sympy cannot use math module functions, but it does have its own equivalents. Sympy supports
        logic operators ~ (not), & (and), | (or), and ^ (xor), so the equation can be multiple equations
        separated by operators.

        If a callable is given to parameter zones, it must take a MesaData object and return a boolean
        mask for what zones to output in a numpy array.


        Example
        --------
        plotter = AutoMESAPlotter( 
                       AutoMESA_run_dir="path/to/AutoMESA_run", 
                       AutoMESA_gridfile="path/to/model_grid_file" )
        data = plotter.get_profile_data_withAutoMESA(
                       gridVarDict={"initial_mass": 1.0, "initial_metallicity": 0.02}, 
                       profile="profile1.data", 
                       variables=["T", "P"], 
                       zones=None )
        """

        # Raise exception if invalid input datatype
        if not isinstance(profile, str):
            raise TypeError(f"profile must be a string, but was a {type(profile)} instead")

        if isinstance(variables, list):
            if not all(isinstance(v, str) for v in variables):
                raise TypeError(f"variables must be str, list of str, or None")
        elif not isinstance(variables, str|type(None)):
            raise TypeError(f"variables must be str, list of str, or None")
        
        if not isinstance(gridVars, dict):
            raise TypeError(f"gridVars must be a dict, but was a {type(gridVars)} instead")
        else:
            for var, values in gridVars.items():
                if not isinstance(var, str):
                    raise TypeError(f"gridVars key must be str, but was a {type(var)} instead")
                
                if isinstance(values, list):
                    if not all(isinstance(i, int|float) for i in values):
                        raise TypeError(f"gridVars items must be int, float, list of floats, or None")
                elif not isinstance(values, int|float|type(None)):
                    raise TypeError(f"gridVars items must be int, float, list of floats, or None")
        
        if isinstance(zones, list):
            if not all(isinstance(i, int|str) for i in zones):
                raise TypeError(f"zones must be an int, list of ints, list of strs, or None")
        elif not isinstance(zones, int|Callable|type(None)):
            raise TypeError(f"zones must be an int, list of ints, list of strs, or None")

        
        
        # Overwrite parameters with deep copies to avoid changing the user's variables.
        gridVars = copy.deepcopy(gridVars)
        zones = copy.deepcopy(zones)
        variables = copy.deepcopy(variables)


        # Put int and str values into lists to avoid continuously checking datatype.
        if isinstance(variables, str): variables = [variables]
        if isinstance(zones, int): zones = [zones]

        use_zone_equation = False
        use_zone_function = False
        if isinstance(zones, list) and all(isinstance(i, str) for i in zones):
            sympy_vars = {}
            for var in zones[1:]:
                sympy_vars[var] = Symbol(var)
            zone_equation = parse_expr(zones[1], local_dict=sympy_vars)
            use_zone_equation = True
        elif isinstance(zones, Callable):
            use_zone_function = True
            
    
        # Check if names in gridVars are valid grid variable names, nicknames, or tied variables.
        # If they are nicknames or tied vars, change the name and values to the corresponding
        # grid variable names and values.
        gridVarGivenNames = list(gridVars.keys())
        gridVarNames = []
        for key in gridVarGivenNames:
            gridVarName = key
            if key not in self.modelGrid.varNames:
                found = False

                # Check if nickname
                for nicknamed_key, nickname in self.modelGrid.nicknamedVars.items():
                    if key == nickname:
                        found = True
                        gridVarName = nicknamed_key
                        gridVars[nicknamed_key] = gridVars.pop(key)
                        break
                
                # Check if tied variable
                if not found:
                    for tiedVar in self.modelGrid.tiedVars:
                        tiedVarName = tiedVar['tied_name']
                        if key == tiedVarName:
                            found = True
                            tiedVarGivenValues = gridVars.pop(key) # Values given by the useer
                            gridVarName = tiedVar['grid_var']
                            varNameIndex = self.modelGrid.varNames.index(gridVarName)
                            if tiedVarGivenValues is None:
                                gridVars[gridVarName] = self.modelGrid.varValues[varNameIndex].tolist()
                                break
                            elif isinstance(tiedVarGivenValues, float):
                                tiedVarValues = [tiedVarValues]
                            
                            gridVarValues = self.modelGrid.varValues[varNameIndex].tolist() # All grid variable values
                            tiedVarValues = [float(tiedVar['rule'].subs(tiedVar['parameter'], value)) for value in gridVarValues] # All tied variable values
                            gridVarGivenValues = [] # Grid variable values that correspond to the values given by the user
                            for value in tiedVarGivenValues:
                                try:
                                    i = tiedVarValues.index(value)
                                except ValueError:
                                    raise ValueError(f"Tied variable {tiedVarName} value {value} does not correspond to any value under the grid variable {gridVarName}.")
                                gridVarGivenValues.append(gridVarValues[i])
                            
                            gridVars[gridVarName] = gridVarGivenValues
                            break
                
                if not found:    
                    raise ValueError(f"gridVars key {key} is invalid")
                
            gridVarNames.append(gridVarName)
            

        # If any gridVars value is None, give it the entire list of values.
        # If a gridVars value is a single int or float, turn it into a list.
        # Put gridVars names and values into lists.
        gridVarValuesList = []
        for gridVarName in gridVarNames:
            if gridVars[gridVarName] is None:
                i = self.modelGrid.varNames.index(gridVarName)
                gridVars[gridVarName] = self.modelGrid.varValues[i].tolist()
            elif isinstance(gridVars[gridVarName], float|int):
                gridVars[gridVarName] = [gridVars[gridVarName]]
            
            gridVarValuesList.append(gridVars[gridVarName])
        

        # Iterate through every combination of each value given in gridVars, which is every model specified by gridVars.
        # For each model, get the values for each variable at the specified zones and put them in dictionaries.
        # Put these dictionaries together in the complete dictionary under the model variables and values
        mult_models_mult_vars = {}
        check_variables = True
        for gridVarValues in iter_product(*gridVarValuesList):
            model = self.get_model(dict(zip(gridVarNames, gridVarValues)))
            mr_profile = mr.MesaData(f'{self.localAutoMESADir}{self.runName}{model.modelName}/LOGS/{profile}')

            # In the first iteration, check if variable names are found in profile.
            # If variables is None, 
            if check_variables:
                if variables is None:
                    variables = list(mr_profile.bulk_data.keys())
                else:
                    for variable in variables:
                        if not mr_profile.in_data(variable):
                            raise ValueError(f"Variable {variable} not found in profile")
                if use_zone_equation:
                    for var in sympy_vars:
                        if not mr_profile.in_data(var):
                            raise ValueError(f'Variable in zones tuple not found in profile')
                check_variables = False
            
            use_zone_mask = True
            if use_zone_equation:
                zone_mask = np.zeroes(len(mr_profile.data('zone')), dtype=bool)
                for i in range(len(mr_profile.data('zone'))):
                    sympy_var_values = {}
                    for var in sympy_vars:
                        sympy_var_values[var] = mr_profile.data(var)
                    zone_mask[i] = zone_equation.subs(sympy_var_values)
            elif use_zone_function:
                zone_mask = zones(mr_profile)
            else:
                use_zone_mask = False

            values_dicts = {}
            for variable in variables:
                # In the first iteration, check if variable names are valid.
            
                if zones is None:
                    values = mr_profile.data(variable)
                elif use_zone_mask:
                    values = mr_profile.data(variable)[zone_mask]
                else:
                    values = np.empty(len(zones))
                    for i in range(len(zones)):
                        values[i] = mr_profile.data(variable)[zones[i]]
                    
                values_dicts[variable] = values
            
            mult_models_mult_vars[tuple(zip(gridVarGivenNames, gridVarValues))] = values_dicts
        


        # Returns un-simplified data if specified by parameter.
        if not simplify:
            return mult_models_mult_vars

        # Simplify return value if possible to avoid unneccessary nesting.
        if len(mult_models_mult_vars) == 1:
            one_model_mult_vars = dict(mult_models_mult_vars.popitem())
            if len(variables) == 1:
                one_model_one_var = one_model_mult_vars.popitem()[1]
                if len(zones) == 1:
                    one_value = one_model_one_var[0]
                    return one_value
                else:
                    return one_model_one_var
            else:
                return one_model_mult_vars
        elif len(variables) == 1:
            mult_models_one_var = {}
            for key, value in mult_models_mult_vars.items():
                mult_models_one_var[key] = value.popitem()[1]
            return mult_models_one_var
        else:
            return mult_models_mult_vars
                
        ...