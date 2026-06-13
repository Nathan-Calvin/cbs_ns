# Author:  Jess Vriesema, Nathan Steenwyk and Ethan Webber
# Date:    4 June 2026
# Purpose: Define plotting tools for use with AutoMESA.


import numpy as np
import matplotlib.pyplot as plt
import mesa_reader as mr
from itertools import product
from MESA_models import MESAModelGrid
from MESA_models import MESAModel
from dataclasses import dataclass



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
                 indices: int|list[int]=None|None) -> np.ndarray:


Examples of loading data from MESA output files:
-----------------------------------------------
    # Grabs temp and pressure data from profile1.data for the model with 
    #   mdot=1.0 and Ltrans=0.02 in the model grid:
    data = plotter.get_profile_data_with_auto_mesa( 
                    gridVars=(("mdot", 1.0), ("Ltrans", 0.02)),
                    profile="profile_600.data",
                    variables=["T", "P"],
                    indices=None )
    # Plot P vs T for that model:
    plt.plot( data["T"], data["P"] )

    # Grabs the surface temp and pressure for multiple models:
    data = plotter.get_profile_data_with_auto_mesa(
                    gridVars=(("mdot", None), ("Ltrans", 0.02)),
                    profile="profile_600.data",
                    variables=["T", "P"],
                    indices=None )
    # Plot Psurf vs Tsurf for each model on the same plot for different mdot values:
    for mdot in data.keys():
        plt.plot( data[mdot]["T"], data[mdot]["P"] )

    # Grabs the star's outer radius at the last profile for a whole 2D grid of models:
    data = plotter.get_profile_data_with_auto_mesa(
                    gridVars=(("mdot", None), ("Ltrans", None)),
                    profile="profile_end.data",
                    variables=["R"],
                    indices=1 )
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
                raise ValueError("Value not found under variable name.")
            
            modelIndex.append(index)

        return self.modelGrid.models[tuple(modelIndex)]




    # For Ethan
    def get_data_simple(self, 
                 autoMESADir:str, 
                 modelName:str, 
                 profile:str, 
                 variables: str|list[str]=None, 
                 indices: int|list[int]|None=None) -> np.ndarray:
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
        
        indices (int, list of ints, or None): 
           The indices of the data series to return. If indices is None, 
           then return everything. If it is just an integer, return just 
           one value. If indices is a list of integers, then return only the 
           data at those indices. All indices should be non-negative integers. 
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
        if not isinstance(indices, (int, list, type(None))):
            raise TypeError(f"indices must be an int, a list of ints, or None, but was a {type(indices)} instead.")

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
                    #If indices is not None, then return only the data at those specified indices.
                    if isinstance(indices, list) or indices is None:
                        # We need to return a list of values at those indices
                        extracted_data[var] = data_array[indices] if indices is not None else data_array
                    elif isinstance(indices, int):
                        # We need to return just one value at that index
                         extracted_data[var] = data_array[indices] 
                    else:
                        # raise an error! indices should be an int, a list of ints, or None.
                        raise TypeError(f"indices must be an int, a list of ints, or None, but was a {type(indices)} instead.")
    
            #Returns a dictionary that maps each variable to its data. 
            return extracted_data 
        elif isinstance(variables, str):
             #Checks if the single string exsits in the MESA file output file. 
            if not md.in_data(variables):
                raise ValueError(f"Variable '{variables}' not found in MESA output file '{profile}' for model '{modelName}' in AutoMESA directory '{autoMESADir}'.")
            else:
                data_array = md.data(variables)
                #If indices is not None, then return only the data at those specified indices.
                if indices is not None: 
                    data_array = data_array[indices]
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
                 variables: str|list[str]=None, 
                 indices: int|list[int]|None=None) -> np.ndarray:
        """Returns the data series at all levels for a given variable from 
        a given MESA output file. 

        PARAMETERS

        gridVars dict[str, float|list[float]|None]:
            A dictionary mapping each grid variable to its value for the 
            model that we want to get the data for. For example, if the 
            grid variables are "initial_mass" and "initial_metallicity", 
            then this dictionary might look like {"initial_mass": 1.0, 
            "initial_metallicity": 0.02}. If there is no model in the model 
            grid with those grid variable values, then this method should 
            raise an error.
            TODO: We should allow for the grid variables to be specified
            by their nicknames, if they exist. 
            TODO: We should also allow for the grid variables to be specified by
            their TiedVariable names, if they exist. 
            TODO: We could also allow gridVarDict to be a list of dictionaries
            or a dictionary mapping each grid variable to a list of values.
            In this case, we would return a list of data series, one for 
            each model specified by the dictionaries. But maybe we should 
            just have the user call this method multiple times if they want 
            data for multiple models, since that would be simpler to implement 
            and use.

        variables (str or list of str): 
            The variable(s) to get the data for. 
            If variables is a string, then return the data series for that variable.
            If variables is a list of strings, then return a dictionary mapping 
            each variable to its data.
        
        indices (int, list of ints, or None): 
           The indices of the data series to return. If indices is None, 
           then return everything. If it is just an integer, return just 
           one value. If indices is a list of integers, then return only the 
           data at those indices. All indices should be non-negative integers. 
           If any index is out of bounds for the data series, then raise 
           an error.
        
        RETURNS
        float, if variables is a string and indices is an integer and the 
            data series has only one value.
        np.ndarray[float], if only one model is specified and variables is a string.
        dict[str, np.ndarray[float]], if only one model is specified and variables 
            is a list of strings. The keys of the dictionary are the variables and the 
            values are the data series for those variables.
        dict[tuple[tuple[str, float]], np.ndarray[float]], if multiple 
            models are specified and only one variable is requested. 
            The keys of the dictionary are a tuple of tuples representing 
            the model identifiers. For example, if the grid variables are 
            "initial_mass" and "initial_metallicity", then the keys might 
            look like (("initial_mass", 1.0), ("initial_metallicity", 0.02)). 
            The values of the dictionary are the data series for that variable for each model.
        dict[tuple[tuple[str, float]], dict[str, np.ndarray[float]]], if multiple 
            models are specified and variables is a list of strings. 
            The keys of the outer (first) dictionary are the model identifiers 
            (e.g. a string representation of the grid variables and their values) 
            and the inner (second) dictionary mapping each variable to its data 
            series for that model.
        """

        # Raise exception if invalid input
        if not isinstance(profile, str):
            raise TypeError(f"profile must be a string, but was a {type(profile)} instead")
        if not isinstance(indices, (int, list, type(None))):
            raise TypeError(f"indices must be an int, a list of ints, or None, but was a {type(indices)} instead")

        if isinstance(variables, list):
            if not all(isinstance(v, str) for v in variables):
                raise TypeError("variables must be str or list of str")
        elif not isinstance(variables, str):
            raise TypeError("variables must be str or list of str")
        
        if not isinstance(gridVars, dict):
            raise TypeError("gridVars must be a dict")
        else:
            for var, values in gridVars.items():
                if not isinstance(var, str):
                    raise TypeError("gridVars key must be str")
                
                if isinstance(values, list):
                    if not all(isinstance(i, float) for i in values):
                        raise TypeError("gridVars items must be float, list of floats, or None")
                elif not isinstance(values, float|type(None)):
                    raise TypeError("gridVars items must be float, list of floats, or None")

        # Put indices and variables into lists to avoid needing to continuously check datatype.
        if isinstance(indices,int): indices = [indices]
        if isinstance(variables,str): variables = [variables]
    
        # Check if names in gridVars are valid. If any are not, check if they are nicknames and change them to valid names.
        # If they are invalid, and aren't nicknames, raise an exception.
        keys = list(gridVars.keys())
        for key in keys:
            if key not in self.modelGrid.varNames:
                found = False
                for nicknamed_key, nickname in self.modelGrid.nicknamedVars.items():
                    if key == nickname:
                        gridVars[nicknamed_key] = gridVars.pop(key)
                        found = True
                        break
                if not found:
                    raise TypeError(f"gridVars name {key} is invalid")
            
        # Purpose 1) If any gridVars value is None, give it the entire list of values.
        # If gridVars value is a single float, not a list of floats, turn it into a list.
        # Purpose 2) Put gridVars keys and values into lists.
        gridVarNames = list(gridVars.keys())
        gridVarValuesList = []
        for gridVarName in gridVarNames:
            if gridVars[gridVarName] is None:
                i = self.modelGrid.varNames.index(gridVarName)
                gridVars[gridVarName] = self.modelGrid.varValues[i].tolist()
            elif isinstance(gridVars[gridVarName], float):
                gridVars[gridVarName] = [gridVarName]
            
            gridVarValuesList.append(gridVars[gridVarName])
        
        # Iterate through every combination of each value given in gridVars, which is every model specified by gridVars.
        # For each model, get the values for each variable at the specified indices and put them in dictionaries.
        # Put these dictionaries together in the complete dictionary under the model variables and values
        mult_models_mult_vars = {}
        for gridVarValues in product(*gridVarValuesList):
            model = self.get_model(dict(zip(gridVarNames, gridVarValues)))
            mr_profile = mr.MesaData(f'{self.localAutoMESADir}{self.runName}{model.modelName}/LOGS/{profile}')
            values_dicts = {}
            for variable in variables:
                if indices is None:
                    values = mr_profile.data(variable)
                else:
                    values = np.empty(len(indices))
                    for i in range(len(indices)):
                        values[i] = mr_profile.data(variable)[indices[i]]
                values_dicts[variable] = values

            for i in range(len(gridVarNames)):
                nickname = self.modelGrid.nicknamedVars.get(gridVarNames[i])
                if nickname is not None:
                    gridVarNames[i] = nickname
            
            mult_models_mult_vars[tuple(zip(gridVarNames, gridVarValues))] = values_dicts
        
        if len(mult_models_mult_vars) == 1:
            one_model_mult_vars = dict(mult_models_mult_vars.popitem())
            if len(variables) == 1:
                one_model_one_var = one_model_mult_vars.popitem()[1]
                if len(indices) == 1:
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
                
        # Example of how this method could be used:
        # plotter = AutoMESAPlotter( 
        #               AutoMESA_run_dir="path/to/AutoMESA_run", 
        #               AutoMESA_gridfile="path/to/model_grid_file" )
        # data = plotter.get_profile_data_withAutoMESA(
        #               gridVarDict={"initial_mass": 1.0, "initial_metallicity": 0.02}, 
        #               profile="profile1.data", 
        #               variables=["T", "P"], 
        #               indices=None )
        ...