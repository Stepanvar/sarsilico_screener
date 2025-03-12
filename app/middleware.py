# middleware.py
from django.http import JsonResponse
from django.shortcuts import redirect
from django.urls import resolve


class PipelineValidationMiddleware:
    """Ensures screening pipeline sequence is followed"""
    
    def __init__(self, get_response):
        self.get_response = get_response
        self.pipeline_steps = [
            'target_selection',
            'input',
            'validation',
            'similarity_check',
            'docking_progress',
            'results'
        ]

    def __call__(self, request):
        response = self.get_response(request)
        return response

    def process_view(self, request, view_func, view_args, view_kwargs):
        current_step = resolve(request.path_info).url_name
        
        if current_step in self.pipeline_steps:
            if not self.validate_pipeline_progression(request, current_step):
                return redirect('target_selection')

    def validate_pipeline_progression(self, request, current_step):
        """Check session for valid pipeline state"""
        session_steps = request.session.get('pipeline_steps', [])
        
        # Allow returning to previous steps
        if current_step in session_steps:
            return True
            
        # Check sequential progression
        expected_step = self.pipeline_steps[len(session_steps)]
        if current_step != expected_step:
            return False
            
        request.session['pipeline_steps'] = session_steps + [current_step]
        return True

class SecurityHeadersMiddleware:
    """Adds security headers for molecular data"""
    
    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        response = self.get_response(request)
        response['Content-Security-Policy'] = "default-src 'self'"
        response['X-Content-Type-Options'] = 'nosniff'
        return response