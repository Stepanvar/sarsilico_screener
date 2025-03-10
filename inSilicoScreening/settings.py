import os
from pathlib import Path
import posixpath
from decouple import config
from sqlalchemy import false
from celery.schedules import crontab

BASE_DIR = Path(__file__).resolve().parent.parent

# Secret Key (Ensure to keep this secret in production)
SECRET_KEY = '4f0a8e4e-007c-4584-8913-2967f102f503'
INSTALLED_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django_bootstrap5',
    'app',  # Main application
]

# settings.py
ALLOWED_HOSTS = [
    '10.2.0.32',   # Server IP
    '127.0.0.1',    # Localhost IPv4
    'localhost',     # Localhost alias
    'zenome-school-2' # Server hostname (if applicable)
]


MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
]

ROOT_URLCONF = 'inSilicoScreening.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [BASE_DIR / 'app' / 'templates'],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
                'django.template.context_processors.static',
            ],
        },
    },
]

WSGI_APPLICATION = 'inSilicoScreening.wsgi.application'

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': BASE_DIR / 'db.sqlite3',
    }
}

AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]

# Internationalization
LANGUAGE_CODE = 'en-us'
# Static and Media Files
STATIC_URL = '/static/'
STATICFILES_DIRS = [BASE_DIR / 'app' / 'static']
STATIC_ROOT = posixpath.join(*(BASE_DIR.parts + ('static',)))

MEDIA_URL = '/media/'
MEDIA_ROOT = BASE_DIR / 'media'

DEFAULT_AUTO_FIELD = 'django.db.models.BigAutoField'

CELERY_BROKER_URL = 'redis://localhost:6379/0'
CELERY_RESULT_BACKEND = 'redis://localhost:6379/0'
CELERY_ACCEPT_CONTENT = ['json']
CELERY_TASK_SERIALIZER = 'json'
CELERY_RESULT_SERIALIZER = 'json'
CSRF_TRUSTED_ORIGINS = [
    'http://127.0.0.1:8000',
    'http://127.0.0.1:80',
    'http://10.2.0.32:8000',
    'http://10.2.0.32:8001',
    'http://127.0.0.1:8001',
]

# Redirect URLs after login/logout
LOGIN_REDIRECT_URL = 'home'
LOGOUT_REDIRECT_URL = 'home'

# settings.py
FILE_UPLOAD_MAX_MEMORY_SIZE = 26214400  # 25MB
DATA_UPLOAD_MAX_MEMORY_SIZE = 26214400   # 25MB


# Security Settings
SECURE_HSTS_SECONDS = 3600  # Start with 1 hour, then increase
SECURE_HSTS_INCLUDE_SUBDOMAINS = True
SECURE_HSTS_PRELOAD = True
SECURE_SSL_REDIRECT = False  # Redirect HTTP to HTTPS
SESSION_COOKIE_SECURE = False
CSRF_COOKIE_SECURE = False
DEBUG = True  # Never run with DEBUG=True in production!

EMAIL_BACKEND = 'django.core.mail.backends.console.EmailBackend'

# settings.py
LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
        },
    },
    'root': {
        'handlers': ['console'],
        'level': 'DEBUG',
    },
}
