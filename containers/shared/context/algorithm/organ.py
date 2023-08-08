import dataclasses
import json
import logging
import typing as t
from pathlib import Path

Organ = t.TypeVar("Organ")


@dataclasses.dataclass
class OrganLookup(t.Generic[Organ]):
    mapping_file: Path

    def get(self, id: str) -> Organ:
        for key, organ in self.__get_options():
            if key.lower() == id.lower():
                return organ
        raise ValueError(f"Organ '{id}' is not supported")

    def get_builtin_options(self) -> t.Iterable[t.Tuple[str, Organ]]:
        return []

    def from_raw(self, raw: t.Any) -> Organ:
        return raw

    def __get_options(self) -> t.Iterable[t.Tuple[str, Organ]]:
        yield from self.get_builtin_options()
        try:
            for key, value in self.__load_mapping_file():
                yield key, self.from_raw(value)
        except ValueError:
            logging.warn(f"Invalid format of organ mapping file '{self.mapping_file}'")

    def __load_mapping_file(self) -> t.Iterable[t.Tuple[str, t.Any]]:
        if not self.mapping_file.exists() or not self.mapping_file.is_file():
            return
        with open(self.mapping_file) as file:
            data = json.load(file)
        yield from data.items() if isinstance(data, dict) else data
