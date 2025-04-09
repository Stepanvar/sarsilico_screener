# settings.py 
import os
from pathlib import Path
import os
from pathlib import Path
import posixpath

# Base directory
BASE_DIR = Path(__file__).resolve().parent.parent

# Debug mode (set to False in production)
DEBUG = True

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

# Database (using SQLite for simplicity)
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': BASE_DIR / 'db.sqlite3',
    }
}


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

WSGI_APPLICATION = 'inSilicoScreening.wsgi.application'

# Static and Media Files
STATIC_URL = '/static/'  # URL prefix for static files
STATIC_ROOT = BASE_DIR / 'staticfiles'  # Collected files directory
STATICFILES_DIRS = [
    BASE_DIR / 'static',  # Your main static files directory
    BASE_DIR / 'app/static/app',  # App-specific static files (if needed)
]

MEDIA_URL = '/media/'
MEDIA_ROOT = BASE_DIR / 'media'

DEFAULT_AUTO_FIELD = 'django.db.models.BigAutoField'

CELERY_BROKER_URL = 'redis://localhost:6379/0'
CELERY_RESULT_BACKEND = 'redis://localhost:6379/0'
CELERY_ACCEPT_CONTENT = ['json']
CELERY_TASK_SERIALIZER = 'json'
CELERY_RESULT_SERIALIZER = 'json'
CORS_ALLOWED_ORIGINS = [
]
CSRF_TRUSTED_ORIGINS = [
    'http://10.2.0.32:8001',
    'http://127.0.0.1:8000',
    'http://127.0.0.1:80',
    'http://10.2.0.32:8000',
    'http://10.2.0.32:8001',
    'http://127.0.0.1:8001',
    'https://*.chemicalize.com',  # For MarvinJS web services
    'https://marvinjs.chemicalize.com'
]

# Redirect URLs after login/logout
LOGIN_REDIRECT_URL = 'index'
LOGOUT_REDIRECT_URL = 'index'

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

AUTH_PASSWORD_VALIDATORS = [
    {'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator'},
    {'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator'},
    {'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator'},
    {'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator'},
]

LANGUAGE_CODE = 'en-us'
TIME_ZONE = 'UTC'
USE_I18N = True
USE_TZ = True

# Celery Configuration
CELERY_CACHE_BACKEND = 'django-cache'

# Custom project settings
PDB_VALIDATION_RESOLUTION = 3.0
MAX_SMILES_LENGTH = 500
MAX_SEQUENCE_LENGTH = 2000
# inSilicoScreening/settings.py
COVIDRUG_API_ENDPOINT = "https://covirus.cc/api/v2"
# Add these to Django settings
SUBPROCESS_TIMEOUT = 300  # 5 minutes
VINA_CONFIG = {
    'BOX_SIZE': [20, 20, 20],
    'EXHAUSTIVENESS': 32,
    'CPU_CORES': 4,
    'RANDOM_SEED': 42
}
ALPHAFOLD_TIMEOUT = 600  # 10 minutes
MAX_PDB_RESOLUTION = 3.0
# settings.py
TTD_API_ENDPOINT = "https://api.covirus.cc/drugs/v1"
TTD_API_KEY = os.getenv("TTD_API_KEY")
JSON_ENCODER = 'app.utils.ChemistryJSONEncoder'

"""
Project configuration (insilico.txt dependencies)
"""

# API Endpoints
TTD_API_URL = "https://api.pharmadb.org/ttd/v1" 
COVIDRUG_SCRAPE_URL = "https://covdrug.org/drug"
ALPHAFOLD_DATA_DIR = "/home/s-zuev/alphafold"
# File Paths
ALPHAFOLD_DOCKER_IMAGE = "alphafold:latest"
VINA_BOX_SIZE = [20, 20, 20]  # Docking grid dimensions
VINA_EXHAUSTIVENESS = 8       # Search intensity
VINA_CPU_CORES = 4            # Processing threads
VINA_NUM_POSES = 10           # Output poses
# Validation Thresholds
SIMILARITY_THRESHOLD = 0.85  # Section 3
MAX_PDB_RESOLUTION = 3.0     # Section 2

# settings.py
MIGRATION_MODULES = {
    'app': 'app.migrations'  # Explicitly set migration path
}
AUTH_USER_MODEL = 'auth.User'  # Explicitly define user model

# Other settings as needed (e.g., ALPHAFOLD_DOCKER_IMAGE, SUBPROCESS_TIMEOUT, etc.)
# settings.py

# CSP_FONT_SRC = ("'self'", "https://cdnjs.cloudflare.com")
# CSP_IMG_SRC = ("'self'", "data:")
# CSP_CONNECT_SRC = ("'self'", "https://covirus.cc", "https://api.pharmadb.org")

# CSP_DEFAULT_SRC = ("'self'", "'unsafe-inline'", "'unsafe-eval'",    # Needed for JSME's code generation
#     "'strict-dynamic'")
# # settings.py
# CSP_STYLE_SRC_ELEM = ("'self'", "https://cdn.jsdelivr.net")
# CSP_SCRIPT_SRC = [
#     "'self'",
#     "https://cdn.jsdelivr.net",
#     "https://yastatic.net",  # Add Yandex CDN
#     "'strict-dynamic'",
#     "'nonce-{{ request.csp_nonce }}'"  # Use Django's nonce
# ]

# CSP_STYLE_SRC = [
#     "'self'",
#     "https://cdn.jsdelivr.net",
#     "https://cdnjs.cloudflare.com",
#     "'unsafe-inline'"  # Required for Bootstrap/jQuery inline styles
# ]

# CSP_INCLUDE_NONCE_IN = ['script-src']
# # settings.py
STATICFILES_STORAGE = 'django.contrib.staticfiles.storage.ManifestStaticFilesStorage'
