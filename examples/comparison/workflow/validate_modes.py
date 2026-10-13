"""Check four mock mode declarations."""
from pathlib import Path
import json

from base_edit_safety.population.common import PopulationError
from base_edit_safety.population.comparison_provenance import resolve_study
from base_edit_safety.population.comparison_questions import prepare_question_set, evaluate_requirements


def main():
    base = Path(__file__).resolve().parent / 'mode-declarations'
    results = []
    for directory in sorted(path for path in base.iterdir() if path.is_dir()):
        expected = json.loads((directory / 'expected.json').read_text())['status']
        try:
            study = resolve_study(directory / 'study.json')
            questions = prepare_question_set(study, [])
            selected = [key for key, value in questions.eligibility.items() if value['eligible']]
            result = evaluate_requirements(questions, selected)
            status = 'complete' if result['unweighted_completion'] else 'ineligible'
            assert all(row['evidence_status'] == 'no_measurements' for row in result['questions'].values())
        except PopulationError:
            status = 'invalid'
        if status != expected:
            raise RuntimeError('Artificial mode expectation did not match')
        results.append({'case': directory.name, 'status': status})
    print(json.dumps({'scope': 'artificial declarations and eligibility only', 'checks': len(results),
                      'results': results}, sort_keys=True))


if __name__ == '__main__':
    main()
