from django.conf import settings

def settings_access(request):
    return {'settings': settings}