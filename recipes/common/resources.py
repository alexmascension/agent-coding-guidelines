"""Locate files shipped inside the recipes package (R scripts), from checkout or wheel."""

from importlib import resources
from pathlib import Path


def package_script(package: str, name: str) -> Path:
    """Return the path of a script distributed with ``package`` (e.g. ``"recipes.scrnaseq"``).

    Raises
    ------
    FileNotFoundError
        If the script is not part of the installed package (packaging error).
    """
    path = Path(str(resources.files(package).joinpath(name)))
    if not path.is_file():
        raise FileNotFoundError(
            f"{name} not found in package {package}; check package data in pyproject.toml"
        )
    return path
