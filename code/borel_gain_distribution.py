"""Distribution of eligible root gains over ALL three-edge-colourings.

For each folded diagram, every proper three-edge-colouring of the closed
graph H (up to global colour names) has two eligible matchings: the two
colour classes not containing the closing edge.  Each has a root gain in
GF(q)^*.  This script records, for stored canonical failures and for
control diagrams, how those gains are distributed, in particular how many
are one.  It is a diagnostic for choosing a proof strategy (counting versus
construction); it proves nothing uniform.
"""

import argparse
import json
import time
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

from borel_exchange_stress import LogField
from borel_folded_experiment import from_point
from borel_gain_lattice import chord_logs, root_row
from borel_kempe_closure import all_three_colourings
from borel_kempe_exchange import closed_graph, normalized
from borel_reversal_experiment import PointDiagram
from borel_three_exchange import audit_matching

FIELDS = {}


def field_for(q):
    if q not in FIELDS:
        FIELDS[q] = LogField(q)
    return FIELDS[q]


def diagram_stats(q, t, c, lift_checks=4):
    field = field_for(q)
    cycle = field.cycle(t)
    assert len(cycle) == q+1
    point = PointDiagram(field, t, c, cycle)
    folded = from_point(point)
    graph, natural = closed_graph(folded)
    logs = chord_logs(folded)
    start = time.time()
    states = all_three_colourings(graph)
    seconds = time.time()-start
    exponents = Counter()
    per_colouring = Counter()   # number of unit gains among the two eligible classes
    rows = {}                   # signed root row -> (exponent, multiplicity)
    natural_gains = None
    natural_state = normalized(natural)
    checked = 0
    for state in sorted(states):
        gains = []
        for excluded in sorted({0, 1, 2}-{state[-1]}):
            raw = root_row(folded, graph, state, excluded)
            exponent = sum(a*b for a, b in zip(raw, logs)) % (q-1)
            gains.append(exponent)
            exponents[exponent] += 1
            key = tuple(raw)
            if key in rows:
                assert rows[key][0] == exponent
                rows[key][1] += 1
            else:
                rows[key] = [exponent, 1]
            if checked < lift_checks:
                selected = frozenset(e for e, col in enumerate(state) if col == excluded)
                gain = field.exp[exponent]
                assert folded.root_path(selected)['gain'] == gain
                certificate, _ = audit_matching(point, folded, selected, permutations=True)
                assert certificate['root_gain'] == gain
                assert certificate['colours'] == (gain != 1)
                checked += 1
        per_colouring[sum(g == 0 for g in gains)] += 1
        if tuple(state) == natural_state:
            natural_gains = gains
    total = sum(exponents.values())
    unit = exponents[0]
    row_unit = sum(1 for e, _ in rows.values() if e == 0)
    multiplicities = sorted((m for _, m in rows.values()), reverse=True)
    row_exponents = Counter(e for e, _ in rows.values())
    distinct = len(exponents)
    largest = exponents.most_common(1)[0]
    return {'q': q, 't': t, 'c': c, 'partner': c ^ t, 'h': folded.h,
            'colourings_mod_global_names': len(states),
            'eligible_matchings': total,
            'unit_gain_matchings': unit,
            'unit_fraction': unit/total,
            'uniform_reference': 1/(q-1),
            'ratio_to_uniform': (unit/total)*(q-1),
            'colourings_with_both_unit': per_colouring[2],
            'colourings_with_one_unit': per_colouring[1],
            'colourings_with_no_unit': per_colouring[0],
            'distinct_root_paths': len(rows),
            'unit_gain_root_paths': row_unit,
            'root_path_unit_fraction': row_unit/len(rows),
            'root_path_ratio_to_uniform': row_unit/len(rows)*(q-1),
            'largest_root_path_multiplicities': multiplicities[:5],
            'root_path_histogram': dict(sorted(row_exponents.items())),
            'distinct_gain_values': distinct,
            'most_common_gain_exponent': largest[0], 'most_common_count': largest[1],
            'natural_colouring_unit_gains': (None if natural_gains is None
                                             else sum(g == 0 for g in natural_gains)),
            'independent_lift_checks': checked,
            'enumeration_seconds': round(seconds, 1),
            'histogram': dict(sorted(exponents.items()))}


def worker(args):
    q, t, c, tag = args
    result = diagram_stats(q, t, c)
    result['tag'] = tag
    print('DONE', tag, q, t, c, 'colourings', result['colourings_mod_global_names'],
          'unit', result['unit_gain_matchings'], '/', result['eligible_matchings'],
          'ratio_to_uniform %.2f' % result['ratio_to_uniform'],
          'both_unit', result['colourings_with_both_unit'],
          'paths', result['distinct_root_paths'], 'unit_paths', result['unit_gain_root_paths'],
          'natural_unit', result['natural_colouring_unit_gains'],
          '%.0fs' % result['enumeration_seconds'], flush=True)
    return result


def stored_failures(path):
    jobs = []
    for block in json.loads(Path(path).read_text()):
        seen = set()
        for f in block['canonical_failures']:
            key = (f['t'], min(f['c'], f['c'] ^ f['t']))
            if key in seen:
                continue
            seen.add(key)
            jobs.append((block['q'], f['t'], f['c'], 'failure'))
    return jobs


def controls(q, t, limit):
    """All folded diagrams (one of c, c+t) for one full-cycle trace, up to limit."""
    field = field_for(q)
    assert len(field.cycle(t)) == q+1
    jobs, seen = [], set()
    for c in range(1, q):
        if c == t:
            continue
        key = min(c, c ^ t)
        if key in seen:
            continue
        seen.add(key)
        jobs.append((q, t, c, 'control'))
        if len(jobs) >= limit:
            break
    return jobs


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fields', type=int, nargs='+', default=[64, 128])
    parser.add_argument('--reflected', default='code/borel_reflected_results.json')
    parser.add_argument('--controls', type=int, default=8,
                        help='control diagrams per field for the first stored trace')
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--json-output', type=Path)
    args = parser.parse_args()
    jobs = [j for j in stored_failures(args.reflected) if j[0] in args.fields]
    for q in args.fields:
        t = next(j[1] for j in jobs if j[0] == q)
        jobs += controls(q, t, args.controls)
    jobs.sort(key=lambda j: (j[0], j[3] != 'failure'))
    print('JOBS', len(jobs), flush=True)
    with Pool(args.workers) as pool:
        results = pool.map(worker, jobs, chunksize=1)
    summary = []
    for q in args.fields:
        for tag in ('failure', 'control'):
            rows = [r for r in results if r['q'] == q and r['tag'] == tag]
            if not rows:
                continue
            total = sum(r['eligible_matchings'] for r in rows)
            unit = sum(r['unit_gain_matchings'] for r in rows)
            summary.append({'q': q, 'tag': tag, 'diagrams': len(rows),
                            'colourings': sum(r['colourings_mod_global_names'] for r in rows),
                            'eligible_matchings': total, 'unit_gain_matchings': unit,
                            'unit_fraction': unit/total, 'ratio_to_uniform': unit/total*(q-1),
                            'both_unit_colourings': sum(r['colourings_with_both_unit'] for r in rows),
                            'distinct_root_paths': sum(r['distinct_root_paths'] for r in rows),
                            'unit_gain_root_paths': sum(r['unit_gain_root_paths'] for r in rows)})
    for s in summary:
        print('SUMMARY', s, flush=True)
    if args.json_output:
        args.json_output.write_text(json.dumps(
            {'scope': 'exhaustive three-colouring gain distributions on stored diagrams; '
                      'a strategy diagnostic, not a theorem',
             'summary': summary, 'results': results}, indent=2)+'\n')
