from django.urls import path
from app import views
from django.contrib.auth import views as auth_views
from django.views.generic import TemplateView

urlpatterns = [
    path('', views.home, name='home'),
    path('about/', TemplateView.as_view(
        template_name='app/about.html'
    ), name='about'),
    path('contact/', TemplateView.as_view(
        template_name='app/contact.html'
    ), name='contact'),
    path('login/', auth_views.LoginView.as_view(
            template_name='app/login.html',
            redirect_authenticated_user=True
        ), name='login'),
    path('logout/', auth_views.LogoutView.as_view(
            template_name='app/logout.html'
        ), name='logout'),
    path('smiles-input/', views.smiles_input, name='smiles_input'),
    path('pdb-upload/', views.pdb_upload, name='pdb_upload'),
    path('docking/', views.docking, name='docking'),
    path('docking/results/<uuid:job_id>/', views.docking_results, name='docking_results'),
    path('similarity/', views.similarity_analysis, name='similarity_analysis'),
    path('similarity/results/<uuid:job_id>/', views.similarity_results, name='similarity_results'),
    path('sequence-conversion/', views.sequence_conversion, name='sequence_conversion'),
    path('select/target/', views.target_selection, name='target_selection'),
    path('register/', views.register, name='register'),
]
