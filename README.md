# robokop-genetics
Tools and service wrappers for building Robokop graphs.

#### Caching
To utilize a redis cache, set the following environment variables to match your own redis cache instance:
```
ROBO_GENETICS_CACHE_HOST=localhost
ROBO_GENETICS_CACHE_PORT=6379
ROBO_GENETICS_CACHE_DB=0
ROBO_GENETICS_CACHE_PASSWORD=yourpassword
```

#### Logging and Temporary Files
robokop-genetics depends on a local directory with write permissions for temporary files and logging.

When used in conjunction with robo-commons or robokop-rags, the default robokop home directory will be used. 

For testing or other purposes, set the following environment variable to specify a valid location.
```
ROBO_GENETICS_HOME=/home/example_directory
```
