
import json
import redis
import logging
from robokop_genetics.util import LoggingUtil
from robokop_genetics.simple_graph_components import SimpleEdge, SimpleNode

class GeneticsCache:

    def __init__(self,
                 redis_host="localhost",
                 redis_port=6379,
                 redis_db=0,
                 redis_password="",
                 log_file_path: str = None,
                 prefix: str = ""):
        self.NORMALIZATION_KEY_PREFIX = f'{prefix}normalize-'
        self.logger = LoggingUtil.init_logging(__name__,
                                               logging.INFO,
                                               logFilePath=log_file_path)
        try:
            if redis_password:
                self.redis = redis.Redis(host=redis_host,
                                         port=int(redis_port),
                                         db=int(redis_db),
                                         password=redis_password)
            else:
                self.redis = redis.Redis(host=redis_host,
                                         port=int(redis_port),
                                         db=int(redis_db))
            self.redis.get('x')
            self.logger.info(f"Cache connected to redis at {redis_host}:{redis_port}/{redis_db}")
        except Exception as e:
            self.redis = None
            self.logger.error(e)
            self.logger.error(f"Failed to connect to redis at {redis_host}:{redis_port}/{redis_db}.")

    def set_normalization(self, node_id: str, normalization: tuple):
        normalization_key = f'{self.NORMALIZATION_KEY_PREFIX}{node_id}'
        self.redis.set(normalization_key, json.dumps(normalization))

    def set_batch_normalization(self, batch_dict: dict):
        pipeline = self.redis.pipeline()
        for node_id, normalization in batch_dict.items():
            normalization_key = f'{self.NORMALIZATION_KEY_PREFIX}{node_id}'
            pipeline.set(normalization_key, json.dumps(normalization))
        pipeline.execute()

    def get_normalization(self, node_id: str):
        normalization_key = f'{self.NORMALIZATION_KEY_PREFIX}{node_id}'
        result = self.redis.get(normalization_key)
        normalization = json.loads(result) if result is not None else None
        return normalization

    def get_batch_normalization(self, node_ids: list):
        pipeline = self.redis.pipeline()
        for node_id in node_ids:
            normalization_key = f'{self.NORMALIZATION_KEY_PREFIX}{node_id}'
            pipeline.get(normalization_key)
        results = pipeline.execute()
        for i, result in enumerate(results):
            if results[i] is not None:
                results[i] = json.loads(results[i])
        return results

    def set_service_results(self, service_key: str, results_dict: dict):
        pipeline = self.redis.pipeline()
        for node_id, results in results_dict.items():
            redis_key = f'{service_key}-{node_id}'
            pipeline.set(redis_key, self.__encode_service_results(results))
        pipeline.execute()

    def __encode_service_results(self, service_results: list):
        encoded_results = []
        for (edge, node) in service_results:
            json_node = {"id": node.id, "type": node.type, "name": node.name}
            json_edge = {"source_id": edge.source_id,
                         "target_id": edge.target_id,
                         "provided_by": edge.provided_by,
                         "input_id": edge.input_id,
                         "predicate_id": edge.predicate_id,
                         "predicate_label": edge.predicate_label,
                         "ctime": edge.ctime,
                         "publications": edge.publications,
                         "properties": edge.properties}
            encoded_result = {"edge": json_edge, "node": json_node}
            encoded_results.append(encoded_result)
        return json.dumps(encoded_results)

    def get_service_results(self, service_key: str, node_ids: list):
        pipeline = self.redis.pipeline()
        for node_id in node_ids:
            pipeline.get(f'{service_key}-{node_id}')
            redis_results = pipeline.execute()

        decoded_results = []
        for results in redis_results:
            decoded_results.append(self.__decode_service_results(results) if results is not None else None)
        return decoded_results

    def __decode_service_results(self, redis_results):
        decoded_results = []
        json_object = json.loads(redis_results)
        for result in json_object:
            edge_json = result["edge"]
            node_json = result["node"]
            decoded_results.append((self.__decode_service_edge(edge_json),
                                   self.__decode_service_node(node_json)))
        return decoded_results

    def __decode_service_edge(self, edge_json):
        publications = []
        for p in edge_json['publications']:
            publications.append(p)
        properties = {}
        for key, value in edge_json['properties'].items():
            if key == "distance":
                properties[key] = int(value)
        return SimpleEdge(source_id=edge_json['source_id'],
                          target_id=edge_json['target_id'],
                          provided_by=edge_json['provided_by'],
                          input_id=edge_json['input_id'],
                          predicate_id=edge_json['predicate_id'],
                          predicate_label=edge_json['predicate_label'],
                          ctime=edge_json['ctime'],
                          publications=publications,
                          properties=edge_json['properties'])

    def __decode_service_node(self, node_json):
        # note that right now we're not caching properties or synonyms for service nodes,
        # properties aren't used yet, synonyms will come from normalization after the fact
        return SimpleNode(id=node_json["id"],
                          type=node_json["type"],
                          name=node_json["name"])

    def delete_all_keys_with_prefix(self, prefix: str):
        keys = self.redis.keys(f'{prefix}*')
        if keys:
            self.redis.delete(*keys)
