from dataclasses import dataclass, field
from robokop_genetics.util import Text


@dataclass
class SimpleNode:
    id: str
    type: str
    name: str
    properties: dict = field(default_factory=dict)
    synonyms: set = field(default_factory=set)

    def get_synonyms_by_prefix(self, prefix):
        """Returns curies for any synonym with the input prefix"""
        return set(filter(lambda x: Text.get_curie(x).upper() == prefix.upper(), self.synonyms))


@dataclass
class SimpleEdge:
    source_id: str
    target_id: str
    provided_by: str
    input_id: str
    predicate_id: str
    predicate_label: str
    ctime: int
    publications: list = field(default_factory=list)
    properties: dict = field(default_factory=dict)
