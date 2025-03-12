class SimilarityChecker {
    constructor() {
        this.threshold = 0.85;
    }

    async check(smiles) {
        const response = await fetch('/api/check_similarity', {
            method: 'POST',
            headers: {'Content-Type': 'application/json'},
            body: JSON.stringify({smiles})
        });
        
        return response.json();
    }

    handleResults(data) {
        if(data.similar_drugs.length > 0) {
            this.showWarning(data.similar_drugs);
            return true;
        }
        return false;
    }

    showWarning(drugs) {
        const warningDiv = document.createElement('div');
        warningDiv.className = 'drug-warning';
        warningDiv.innerHTML = `
            <h4>Similar drugs found:</h4>
            <ul>
                ${drugs.map(d => `<li>${d.name} (Similarity: ${d.score})</li>`).join('')}
            </ul>
            <a href="https://covirus.cc/drugs/" target="_blank">
                Check CoviDrug database
            </a>
        `;
        document.body.appendChild(warningDiv);
    }
}
