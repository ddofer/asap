"""
Django settings for cleavepred project.
"""

import sys
import os

### Preparations ###

BASE_DIR = os.path.dirname(__file__)
PROJECT_DIR = os.path.dirname(os.path.dirname(BASE_DIR))
PROJECT_PY_DIR = os.path.join(PROJECT_DIR, 'py')

# In order to load 'cleavepred' module later on
sys.path += [PROJECT_PY_DIR]

### Development vs. Production status ###

DEBUG = True

TEMPLATE_DEBUG = True

ALLOWED_HOSTS = []

### Security ###

SECRET_KEY = 'ij*hpack1#brf--5b_1nd7$cz*8h*y=b!y#_bd48v5kcdg0*vd'

### Localization ###

TIME_ZONE = 'UTC'
LANGUAGE_CODE = 'en-us'

### URLs ###

ROOT_URLCONF = 'urls'

### Context Processors, Middlewares & Apps ###

TEMPLATE_CONTEXT_PROCESSORS = (
    'django.contrib.auth.context_processors.auth',
    'django.core.context_processors.debug',
    'django.core.context_processors.media',
    'django.core.context_processors.static',
    'django.contrib.messages.context_processors.messages',
)

MIDDLEWARE_CLASSES = (
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
)

INSTALLED_APPS = (
    'django.contrib.staticfiles',
    'django.contrib.sessions',
    'django.contrib.contenttypes',
    'django.contrib.auth',
    'django.contrib.admin',
    'django.contrib.humanize',
)

### Templates ###

TEMPLATE_DIRS = (
    os.path.join(BASE_DIR, 'templates'),
)

### Static Files ###

STATIC_URL = '/static/'
