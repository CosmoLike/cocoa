# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

# Jupyter Notebook 7 / jupyter-server configuration
c = get_config()  # get the config object
c.ServerApp.ip = '*'
c.ServerApp.allow_remote_access = True
# do not open a browser window by default when using notebooks
c.ServerApp.open_browser = False
c.ServerApp.root_dir = '/home/whovian/'
# Do not allow Jupyter to be run from the root user inside the container
c.ServerApp.allow_root = False
c.ServerApp.allow_origin = '*'
# ability to change the password at first login time - disabled
c.ServerApp.allow_password_change = False
c.ServerApp.tornado_settings = {
    'headers': {
        'Content-Security-Policy': "frame-ancestors 'self' *"
    }
}
c.ServerApp.max_buffer_size = 70000000000
c.ServerApp.iopub_data_rate_limit = 10000000000
c.ServerApp.shutdown_no_activity_timeout = 0
c.ServerApp.webbrowser_open_new = 2
