import logging
from abc import ABC, abstractmethod
from dataclasses import fields
from typing import TYPE_CHECKING, Any, Mapping, Union

from modular_dales.logging_wrapper import logwrap
from modular_dales.modular.time_dependent_scalars import TimeDependentScalar

logger = logging.getLogger(__name__)
if TYPE_CHECKING:
    import numpy as np
    from modular_dales.Geometry.GridDales import GridDales
    from modular_dales.modular.dales_simulation import dales_simulation
    from modular_dales.vars import VariableDefinition


def set_nml_section(
    nml,
    nml_docs,
    nml_module_setter: str,
    section: str,
    key: str,
    value: Any,
    raise_conflict: bool = False,
) -> None:
    """Set multiple key-value pairs in a namelist section."""
    if section not in nml:
        nml[section] = {}
    if section not in nml_docs:
        nml_docs[section] = {}
    if raise_conflict and key in nml[section] and nml[section][key] != value:
        raise ValueError(
            f"Conflict detected in module '{nml_module_setter}': {key} already exists in section {section} with value {nml[section][key]}, cannot set to {value}! (use manual override)"
        )
    nml[section][key] = value
    nml_docs[section][key] = f"{value} (Set by {nml_module_setter} module)"


class simulation_module(ABC):
    """Base class for DALES simulation modules."""

    def __init__(self, sim: "dales_simulation" = None):
        """Initialize module with reference to parent simulation.

        Args:
            sim: Parent dales_simulation instance providing config, grid, paths, etc.
                 Can be None if module will be added to simulation later.
        """
        self.sim = sim
        self.module_name = None  # Optional name for logging and identification
        self.do_config_done = False
        self.prepare_calculation_done = False
        self.check_settings_done = False
        self.write_files_done = False

    def apply_namelist_from_fields(self) -> None:
        """Apply dataclass field values to the namelist using field metadata.

        Field metadata keys:
            nml: Namelist section (e.g., "NAMSURFACE")
            key: Key in namelist section (defaults to field name)
            required: Whether a value is required (raises if None)
        """

        for dataclass_field in fields(self):
            meta = dataclass_field.metadata or {}
            nml_section = meta.get("nml")
            if not nml_section:
                continue
            key = meta.get("key", dataclass_field.name)
            required = meta.get("required", False)
            if not required:
                logger.info(
                    f"Field '{dataclass_field.name}' in module '{self.module_name}' is optional; skipping if None"
                )
            raise_conflict = meta.get("raise_conflict", False)
            value = getattr(self, dataclass_field.name)

            if isinstance(value, TimeDependentScalar):
                value = value.value_at_start()
            if value is None:
                if required:
                    raise ValueError(
                        f"Module '{self.module_name}' requires '{dataclass_field.name}'"
                    )
                continue
            if isinstance(key, list) or isinstance(nml_section, list):
                if isinstance(value, list):
                    for nml_section, key, value in zip(nml_section, key, value):
                        self.set_nml_section(
                            nml_section,
                            key,
                            value,
                            raise_conflict=raise_conflict,
                        )
                    continue
                elif isinstance(value, dict):
                    for value_key, value_val in value.items():
                        value_section = value_key.split(":")[
                            0
                        ]  # Extract section from "SECTION:KEY"
                        value_key = value_key.split(":")[
                            1
                        ]  # Extract key from "SECTION:KEY"
                        self.set_nml_section(
                            value_section,
                            value_key,
                            value_val,
                            raise_conflict=raise_conflict,
                        )
                else:
                    for nml_section, key in zip(nml_section, key):
                        self.set_nml_section(
                            nml_section,
                            key,
                            value,
                            raise_conflict=raise_conflict,
                        )
                    continue
            self.set_nml_section(
                nml_section,
                key,
                value,
                raise_conflict=raise_conflict,
            )

    def set_nml_section(
        self,
        section: str,
        key: str,
        value: Any,
        raise_conflict: bool = False,
    ) -> None:
        """Helper to set a single namelist section/key from within module methods."""
        set_nml_section(
            self.nml,
            self.nml_docs,
            self.module_name,
            section,
            key,
            value,
            raise_conflict=raise_conflict,
        )

    @property
    def grid(self) -> "GridDales":
        """Access grid from parent simulation."""
        return self.sim.grid if self.sim is not None else None

    @property
    def output_path(self):
        """Access output_path from parent simulation."""
        return self.sim.output_path if self.sim is not None else None

    @property
    def nml(self):
        """Access namelist from parent simulation."""
        return self.sim.nml if self.sim is not None else None

    @property
    def nml_docs(self):
        """Access namelist documentation from parent simulation."""
        return self.sim.nml_docs if self.sim is not None else None

    @property
    def exp_id(self):
        """Access experiment ID from parent simulation."""
        return self.sim.exp_id if self.sim is not None else None

    @property
    def required_folder_list(self):
        """Access required_folder_list from parent simulation.

        This is used by modules that need to register additional
        folders whose contents must be linked or copied to the
        runtime directory (e.g. emissions, radiation tables).
        """

        return self.sim.required_folder_list if self.sim is not None else None

    def retrieve_module(
        self, module_name: Union[str, "simulation_module", type]
    ) -> "simulation_module":
        """Helper to retrieve another module from the parent simulation by name."""
        if self.sim is None:
            raise RuntimeError(
                f"Module '{self.module_name}' is not attached to a simulation!"
            )
        for module in self.sim.modules:
            if module.module_name == module_name:
                return module
            if isinstance(module, module_name):
                return module
        raise ValueError(
            f"Module with name '{module_name}' not found in simulation for retrieval by module '{self.module_name}'"
        )

    def module_exists(self, module_name: Union[str, "simulation_module", type]) -> bool:
        """Check if a module with the given name exists in the parent simulation."""
        if self.sim is None:
            raise RuntimeError(
                f"Module '{self.module_name}' is not attached to a simulation!"
            )
        return any(
            module.module_name == module_name or isinstance(module, module_name)
            for module in self.sim.modules
        )

    @logwrap
    def _initialize_from_sim(self, sim: "dales_simulation"):
        """Initialize module with parent simulation reference."""
        self.sim = sim

    def do_config(self):
        """Configure namelist and settings for this module.

        This is called before prepare_calculation() and should be used by
        configuration modules (like DomainConfigModule, RunConfigModule, etc.)
        to set up namelist parameters.

        Default implementation does nothing. Override in subclasses as needed.
        """
        pass

    def check_settings(self):
        """Check and validate settings for this module.

        This is called after all do_config() methods have completed, allowing
        modules to check for conflicts and validate that the full configuration
        is valid.

        Default implementation does nothing. Override in subclasses as needed.
        """
        pass

    @abstractmethod
    def prepare_calculation(self):
        """Prepare data and setup for calculations."""
        pass

    @abstractmethod
    def write_files(self):
        """Write output files for this module."""
        pass

    def get_timedep_atmosphere_forcings(
        self,
    ) -> Mapping["VariableDefinition", Mapping[float, "np.ndarray"]]:
        """Optional hook: provide time-dependent forcing series for Atmosphere output.

        Expected return format:
            {
                VariableDefinition: {time_seconds: value_or_profile_column, ...},
                ...
            }
        where value_or_profile_column can be a scalar (broadcast over z) or
        a 1D z-profile array.
        :meta private:
        """

        return {}
