# urls.py
from django.urls import path
from . import views
from django.conf import settings
from django.conf.urls.static import static
from django.urls import path

urlpatterns = [
    path('', views.index, name='index'),
    path('login/', views.LoginView.as_view(), name='login'),
    
    path('logout/', views.LogoutView.as_view(), name='logout'),
    
    path('register/', views.CustomRegisterView.as_view(), name='register'),
    path('results/<uuid:job_id>/', views.results_view, name='results'),
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)