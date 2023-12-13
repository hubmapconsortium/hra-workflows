import dataclasses
import json
import logging
import typing as t
from pathlib import Path

Organ = t.TypeVar("Organ")


@dataclasses.dataclass
class OrganLookup(t.Generic[Organ]):
    """Lookup from raw organ name to algorithm specific organ data.

    Attributes:
        mapping_file (Path): Path to file mapping raw organ name to data
    """

    mapping_file: Path

    def get(self, id: str) -> Organ:
        """Get the algorithm specific data for a raw organ name.

        Args:
            id (str): Organ uberon id

        Raises:
            ValueError: If the organ is not supported by the algorithm

        Returns:
            Organ: Algorithm specific data
        """
        for key, organ in self.__get_options():
            if key.lower() == id.lower():
                return organ
        raise ValueError(f"Organ '{id}' is not supported")

    def get_builtin_options(self) -> t.Iterable[t.Tuple[str, Organ]]:
        """Get builtin organ mapping options.

        Returns:
            t.Iterable[t.Tuple[str, Organ]]: Entries mapping organ to data
        """
        return []

    def from_raw(self, raw: t.Any) -> Organ:
        """Convert a raw mapping value to algorithm specific data.
        Can be overridden in subclasses.

        Args:
            raw (t.Any): Raw value from the mapping file

        Returns:
            Organ: Converted organ data
        """
        return raw

    def __get_options(self) -> t.Iterable[t.Tuple[str, Organ]]:
        """Gets all options, builtin and from the mapping file.

        Yields:
            Iterator[t.Tuple[str, Organ]]: Each entry from builtin and the mapping file
        """
        yield from self.get_builtin_options()
        try:
            for key, value in self.__load_mapping_file():
                yield key, self.from_raw(value)
        except ValueError:
            logging.warn(f"Invalid format of organ mapping file '{self.mapping_file}'")

    def __load_mapping_file(self) -> t.Iterable[t.Tuple[str, t.Any]]:
        """Load the mapping json file.

        Yields:
            Iterator[t.Tuple[str, t.Any]]: Each entry in the mapping
        """
        if not self.mapping_file.exists() or not self.mapping_file.is_file():
            return
        with open(self.mapping_file) as file:
            data = json.load(file)
        yield from data.items() if isinstance(data, dict) else data
