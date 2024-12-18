from django.urls import path
from app import views

urlpatterns = [
    path('', views.home_view, name='home'),
    path('about/', views.about_view, name='about'),
    path('contact/', views.contact_view, name='contact'),
    path('login/', views.login_view, name='login'),
    path('logout/', views.logout_view, name='logout'),
    path('input/smiles/', views.smiles_input_view, name='smiles_input'),
    path('input/pdb/', views.pdb_upload_view, name='pdb_upload'),
    path('convert/sequence/', views.sequence_conversion_view, name='sequence_conversion'),
    path('select/target/', views.target_selection_view, name='target_selection'),
    path('docking/results/<str:job_id>/', views.docking_results_view, name='docking_results'),
    path('similarity/', views.similarity_analysis_view, name='similarity_analysis'),
    path('convert_to_smiles/', views.convert_to_smiles, name='convert_to_smiles'),
    # Add more URLs as needed
]
