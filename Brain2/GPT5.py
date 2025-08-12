from brian2 import *
import numpy as np
import matplotlib.pyplot as plt

# Simulation parameters
start_scope()
prefs.codegen.target = 'numpy'
duration = 5 * second  # Training duration
test_duration = 1 * second

# Network size
N_input = 2
N_exc = 100
N_inh = 25

# Input firing rate
rate = 20 * Hz

# Correlation between inputs (0 = uncorrelated, 1 = fully redundant)
input_corr = 0.8  # test with 0.0 and 0.8

# Neuron model
tau = 10*ms
eqs = '''
dv/dt = (-v + I)/tau : 1 (unless refractory)
I : 1
'''

# Input generation with controllable correlation
def generate_correlated_spikes(corr, rate, duration, dt):
    n_steps = int(duration / dt)
    spikes = np.random.rand(N_input, n_steps) < rate * dt
    if corr > 0:
        shared = np.random.rand(n_steps) < rate * dt * corr
        spikes[1] = np.logical_or(spikes[1], shared)
    return spikes

# Create input groups
dt = defaultclock.dt
input_spikes = generate_correlated_spikes(input_corr, rate, duration, dt)
input_times = []
input_indices = []
for i in range(N_input):
    spike_times = np.where(input_spikes[i])[0] * dt
    input_times.extend(spike_times)
    input_indices.extend([i]*len(spike_times))
input = SpikeGeneratorGroup(N_input, input_indices, input_times*second)

# Excitatory and inhibitory neurons
exc = NeuronGroup(N_exc, eqs, threshold='v>1', reset='v=0', refractory=5*ms, method='euler')
inh = NeuronGroup(N_inh, eqs, threshold='v>1', reset='v=0', refractory=5*ms, method='euler')

exc.v = 'rand()'
inh.v = 'rand()'

# Synapses
# Input to Exc
S_input_exc = Synapses(input, exc, 'w : 1', on_pre='I_post += w')
S_input_exc.connect(p=0.5)
S_input_exc.w = '0.2 + 0.1*rand()'

# Input to Inh
S_input_inh = Synapses(input, inh, 'w : 1', on_pre='I_post += w')
S_input_inh.connect(p=0.5)
S_input_inh.w = '0.2 + 0.1*rand()'

# Exc to Exc (recurrent, STDP)
S_ee = Synapses(exc, exc, '''
                w : 1
                dApre/dt = -Apre / tau : 1 (event-driven)
                dApost/dt = -Apost / tau : 1 (event-driven)
                ''',
                on_pre='''
                I_post += w
                Apre += 0.01
                w = clip(w + Apost, 0, 1)
                ''',
                on_post='''
                Apost += 0.01
                w = clip(w + Apre, 0, 1)
                ''')
S_ee.connect(condition='i!=j', p=0.1)
S_ee.w = '0.1 + 0.1*rand()'

# Exc to Inh (feedforward)
S_ei = Synapses(exc, inh, 'w:1', on_pre='I_post += w')
S_ei.connect(p=0.2)
S_ei.w = '0.2 + 0.1*rand()'

# Inh to Exc (feedback inhibition with inhibitory STDP)
S_ie = Synapses(inh, exc, '''
                w : 1
                dApre/dt = -Apre / tau : 1 (event-driven)
                dApost/dt = -Apost / tau : 1 (event-driven)
                ''',
                on_pre='''
                I_post -= w
                Apre += 0.01
                w = clip(w + Apost, 0, 1)
                ''',
                on_post='''
                Apost += 0.01
                w = clip(w + Apre, 0, 1)
                ''')
S_ie.connect(p=0.3)
S_ie.w = '0.2 + 0.1*rand()'

# Monitors
mon_exc_w = StateMonitor(S_input_exc, 'w', record=range(10))  # Monitor 10 input->exc synapses
mon_inh_w = StateMonitor(S_ie, 'w', record=range(10))         # Monitor 10 inh->exc synapses
spike_mon_exc = SpikeMonitor(exc)
spike_mon_inh = SpikeMonitor(inh)

# Run training
run(duration)

# Plot average weights
plt.figure(figsize=(12, 4))
plt.subplot(1, 2, 1)
plt.plot(mon_exc_w.t/ms, np.mean(mon_exc_w.w, axis=0))
plt.title('Avg Excitatory Input Weights')
plt.xlabel('Time (ms)')
plt.ylabel('Weight')

plt.subplot(1, 2, 2)
plt.plot(mon_inh_w.t/ms, np.mean(mon_inh_w.w, axis=0))
plt.title('Avg Inhibitory Weights (Inh -> Exc)')
plt.xlabel('Time (ms)')
plt.ylabel('Weight')
plt.tight_layout()
plt.show()
