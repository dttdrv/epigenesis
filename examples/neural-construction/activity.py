"""Exact, event-driven LIF transmission assay on the constructed contact graph.

Reference equations: NEST iaf_psc_delta, with zero constant current, no lower
voltage bound and refractory input discarded. We implement the equations here;
this does not invoke NEST or establish simulator equivalence. Weights in mV,
delays in ms and neuron parameters are synthetic assay assumptions.
https://nest-simulator.readthedocs.io/en/stable/models/iaf_psc_delta.html
"""
from __future__ import annotations

from collections import defaultdict
import hashlib
import heapq
import json
import math
from pathlib import Path

from brainc._io import read_regular_file


REST_MV = -70.0
THRESHOLD_MV = -55.0
TAU_MS = 10.0
REFRACTORY_MS = 2.0
MAX_EVENTS = 100_000


def simulate(cells, contacts, stimuli, *, duration=20.0, weight=15.0, delay=1.0, max_events=MAX_EVENTS):
    """Numerical kernel; the production assay supplies its graph from a verified trace."""
    if not cells or len(set(cells)) != len(cells) or any(type(cell) is not int or cell < 0 for cell in cells):
        raise ValueError("invalid activity cell identities")
    if any(type(x) not in (int,float) or not math.isfinite(x) for x in (duration,weight,delay)) or duration <= 0 or delay <= 0:
        raise ValueError("invalid activity duration, weight or delay")
    if type(max_events) is not int or not 1 <= max_events <= MAX_EVENTS:
        raise ValueError("invalid activity event budget")
    outgoing = {cell:[] for cell in cells}
    contact_count = 0
    for source,target in contacts:
        if source not in outgoing or target not in outgoing or source == target:
            raise ValueError("activity contact has invalid endpoints")
        outgoing[source].append(target)
        contact_count += 1
        if contact_count > max_events:
            raise ValueError("activity event budget exhausted")
    queue = []
    for time,cell,jump in stimuli:
        if cell not in outgoing or any(type(x) not in (int,float) or not math.isfinite(x) for x in (time,jump)) or not 0 <= time <= duration:
            raise ValueError("invalid external stimulus")
        if len(queue) == max_events:
            raise ValueError("activity event budget exhausted")
        heapq.heappush(queue,(time,cell,jump))
    voltage = dict.fromkeys(cells,REST_MV)
    updated = dict.fromkeys(cells,0.0)
    refractory = dict.fromkeys(cells,0.0)
    spikes, observations = [], []
    processed = 0
    while queue:
        time = queue[0][0]
        inputs = defaultdict(list)
        while queue and queue[0][0] == time:
            _,cell,jump = heapq.heappop(queue)
            inputs[cell].append(jump)
            processed += 1
            if processed > max_events:
                raise ValueError("activity event budget exhausted")
        for cell in sorted(inputs):
            try:
                jump = math.fsum(inputs[cell])
            except OverflowError as failure:
                raise ValueError("nonfinite synaptic input") from failure
            if not math.isfinite(jump):
                raise ValueError("nonfinite synaptic input")
            if time < refractory[cell]:
                observations.append({"time":time,"cell":cell,"jump":jump,"voltage":REST_MV,"refractory":True})
                continue
            voltage[cell] = REST_MV+(voltage[cell]-REST_MV)*math.exp(-(time-updated[cell])/TAU_MS)+jump
            updated[cell] = time
            if not math.isfinite(voltage[cell]):
                raise ValueError("nonfinite membrane potential")
            observations.append({"time":time,"cell":cell,"jump":jump,"voltage":voltage[cell],"refractory":False})
            if voltage[cell] >= THRESHOLD_MV:
                spikes.append([time,cell])
                voltage[cell] = REST_MV
                refractory[cell] = updated[cell] = time+REFRACTORY_MS
                arrival = time+delay
                if arrival <= time:
                    raise ValueError("synaptic delay cannot advance the numerical clock")
                if arrival <= duration:
                    for target in outgoing[cell]:
                        if len(queue)+processed >= max_events:
                            raise ValueError("activity event budget exhausted")
                        heapq.heappush(queue,(arrival,target,weight))
    final = {str(cell):REST_MV+(voltage[cell]-REST_MV)*math.exp(-max(0,duration-updated[cell])/TAU_MS) for cell in sorted(cells)}
    return {"spikes":spikes,"voltages":final,"events":observations}


def assay(directory: Path, verify_trace) -> dict:
    """Probe every developed cell, retaining contact multiplicity and lesion controls."""
    summary = verify_trace(directory)
    raw = read_regular_file(directory/"trace.jsonl",maximum_bytes=64*1024*1024)
    events = [json.loads(line) for line in raw.splitlines()]
    final_step = max(event["step"] for event in events)
    cells = sorted(event["cell"] for event in events if event["event"] == "cell" and event["step"] == final_step)
    contacts = [(event["pre"],event["post"]) for event in events if event["event"] == "connection"]
    trials = []
    passed = not simulate(cells,contacts,[],duration=20)["spikes"]
    for source in cells:
        stimuli = [(1.0,source,THRESHOLD_MV-REST_MV)]
        intact = simulate(cells,contacts,stimuli)
        lesion = simulate(cells,[(pre,post) for pre,post in contacts if pre != source],stimuli)
        blocked = simulate(cells,contacts,stimuli,weight=0.0)
        distance, pending = {source:0}, [source]
        while pending:
            pre = pending.pop(0)
            for edge_pre,post in contacts:
                if edge_pre == pre and post not in distance:
                    distance[post] = distance[pre]+1
                    pending.append(post)
        downstream = sorted({cell for _,cell in intact["spikes"]}-{source})
        expected = {cell for cell,hops in distance.items() if cell != source and 1.0+hops <= 20.0}
        passed &= set(downstream) == expected and lesion["spikes"] == blocked["spikes"] == [[1.0,source]]
        trials.append({"stimulated_cell":source,"intact_spikes":intact["spikes"],"downstream_cells":downstream,
                       "lesioned_spikes":lesion["spikes"],"blocked_spikes":blocked["spikes"]})
    return {"model":"zero-current delta-input LIF", "trace_sha256":hashlib.sha256(raw).hexdigest(),
            "implementation_sha256":hashlib.sha256(read_regular_file(Path(__file__))).hexdigest(),
            "cells":summary["cells"],"contacts":len(contacts),"trials":trials,
            "retained_sites_outside_segments":summary["retained_sites_outside_segments"],
            "assumptions":{"rest_mV":REST_MV,"reset_mV":REST_MV,"threshold_mV":THRESHOLD_MV,
                "tau_ms":TAU_MS,"refractory_ms":REFRACTORY_MS,"contact_jump_mV":15.0,"contact_delay_ms":1.0,
                "duration_ms":20.0,"all_contacts_excitatory":True},
            "transmission_controls_pass":bool(passed),"biological_acceptance":False}


def write_viewer(directory: Path, activity: dict, template: Path, output: Path) -> None:
    events = [json.loads(line) for line in read_regular_file(directory/"trace.jsonl",maximum_bytes=64*1024*1024).splitlines()]
    frames = {}
    for event in events:
        if event["event"] not in ("cell","segment"):
            continue
        frame = frames.setdefault(event["step"],{"time":event["time"],"cells":[],"segments":[]})
        if event["event"] == "cell":
            frame["cells"].append([event["cell"],*event["position"],(3*event["volume"]/(4*math.pi))**(1/3)])
        else:
            frame["segments"].append([event["cell"],event["axon"],*event["proximal"],*event["position"],event["diameter"]])
    data = {"frames":list(frames.values()),"activity":activity}
    text = template.read_text()
    marker = "/* CONSTRUCTION_DATA */ null"
    if text.count(marker) != 1:
        raise ValueError("invalid construction viewer template")
    with output.open("x") as stream:
        stream.write(text.replace(marker,json.dumps(data,allow_nan=False).replace("<","\\u003c")))
