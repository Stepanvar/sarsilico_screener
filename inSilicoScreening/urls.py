# inSilicoScreening/urls.py
from django.contrib import admin
from django.urls import path, include
from django.views.generic import RedirectView
from app import views
from django.contrib.auth import views as auth_views

urlpatterns = [
    # Admin interface
    path('admin/', admin.site.urls),
    
    # Workflow Navigation
    path('', views.screening_view, name='index'),
    path('target_selection/', views.TargetSelectionView.as_view(), name='target_selection'),
    
    # Compound Input Subsystem (grouped under 'input/')
    path('input/', include(([
        path('small-molecule/', views.SmallMoleculeInputView.as_view(), name='small_molecule'),
        path('peptide/', views.PeptideInputView.as_view(), name='peptide'),
        path('biomolecule/', views.BiomoleculeInputView.as_view(), name='biomol'),
    ], 'input'), namespace='input')),
    
    # Validation
    path('validation/<uuid:job_id>/', views.ValidationFeedbackView.as_view(), name='validation_feedback'),
    
    # Workflow Progress Tracking
    path('similarity-check/<uuid:job_id>/', views.SimilarityCheckView.as_view(), name='similarity_check'),
    path('docking/<uuid:job_id>/', views.DockingProgressView.as_view(), name='docking_progress'),
    
    # Results & Visualization (grouped under 'results/')
    path('results/', include(([
        path('<uuid:job_id>/', views.ResultsVisualizationView.as_view(), name='visualization'),
        path('<uuid:job_id>/summary/', views.ResultsSummaryView.as_view(), name='summary'),
    ], 'results'), namespace='results')),
    
    # Authentication endpoints
    path('login/', auth_views.LoginView.as_view(template_name='app/login.html'), name='login'),
    path('logout/', auth_views.LogoutView.as_view(), name='logout'),
    path('register/', views.register, name='register'),
    
    # Optional redirect alias for home
    path('home/', RedirectView.as_view(pattern_name='index', permanent=False)),
]
