# """Tests for SPICE init."""

# import unittest
# import subprocess as sp
# import socket
# import time
# from typing import Tuple, Union
# import os
# import shutil

# import pymongo
# from arrakis.logger import logger

# logger.setLevel("DEBUG")

# from arrakis.init_database import (
#     main,
#     source2beams,
#     ndix_unique,
#     cat2beams,
#     source_database,
#     beam_database,
#     get_catalogue,
#     get_beams,
#     field_database,
# )

# SOURCE_CAT = "data/source_RACS_1237+12A_3165.fits"
# GAUSS = "data/gaussians_RACS_1237+12A_3165.fits"
# EPOCH = "data/arrakis_test"
# FIELD = "data/RACS_1237+12A"

# def start_mongodb(port: int = 27017) -> Tuple[str, int]:
#     """Start a local MongoDB instance."""
#     return
#     if os.path.exists("data/testdb"):
#         shutil.rmtree("data/testdb")

#     os.makedirs("data/testdb")

#     cmd = f"mongod --dbpath data/testdb --port {port} --fork --auth --logpath data/testdb/mongodb.log"
#     logger.debug("Starting mongo...")
#     sp.run(cmd.split(), check=True)
#     logger.debug("Mongo started. Sleeping for 5 seconds...")
#     time.sleep(5)
#     logger.debug("Mongo should be ready.")
#     host = socket.gethostbyname(socket.gethostname())
#     return host, port

# def create_mongo_admin(
#     host: str,
#     port: int=27017,
#     username: str="admin",
#     password: str="admin",
#     ) -> None:

#     """Create an admin user in a local MongoDB instance."""
#     return
#     client = pymongo.MongoClient(host, port)
#     db = client.admin
#     cmd = db.command("createUser", username, pwd=password, roles=["root"])
#     logger.debug(cmd)


# def stop_mongodb() -> None:
#     """Stop a local MongoDB instance."""
#     return
#     cmd = f"mongod --dbpath data/testdb --shutdown"
#     logger.debug("Stopping mongo...")
#     p = sp.Popen(cmd.split(),)
#     logger.debug("Mongo stopped.")

# class TestInit(unittest.TestCase):
#     host, port = start_mongodb()
#     create_mongo_admin(host, port)


#     def test_main(self):
#         main(
#             load=True,
#             islandcat=SOURCE_CAT,
#             compcat=GAUSS,
#             host=self.host,
#             username="admin",
#             password="admin",
#             field=FIELD,
#             epoch=EPOCH,
#             force=True,
#         )
#         assert True


#     logger.debug("BOOP")

#     stop_mongodb()

# if __name__ == '__main__':
#     unittest.main()
