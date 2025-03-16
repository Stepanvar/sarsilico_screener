// Keep JSME applet as global variable
let jsmeApplet = null;

// Initialize when BOTH DOM and JSME are ready
document.addEventListener('DOMContentLoaded', function() {
  // Load JSME dynamically
  const script = document.createElement('script');
  script.src = '/static/jsme/jsme.nocache.js';
  script.onload = function() {
    // JSME-specific initialization (critical!)
    jsmeApplet = new JSApplet.JSME("jsme_container", "1000px", "1000px", {
      "options": "marker"
    });
    
    // Button handler
    document.getElementById('show-smiles').addEventListener('click', function() {
      if(jsmeApplet) {
        alert(jsmeApplet.smiles());
      } else {
        console.error('JSME failed to initialize');
      }
    });
  };
  document.head.appendChild(script);
});
