class ConversionTracker {
    constructor(jobId) {
        this.eventSource = new EventSource(`/api/conversion/${jobId}/status`);
        this._initHandlers();
    }
    
    _initHandlers() {
        this.eventSource.onmessage = (e) => {
            const data = JSON.parse(e.data);
            this._updateProgressBars(data);
            
            if(data.stage === 'complete') {
                window.location.href = '/similarity-check';
            }
        };
    }
}