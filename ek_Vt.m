% === Integrate breathing trace across time to get Vt ===
filtered = highpass(data.cage74_2.beforeCre_pre.voltage, 150);
Vt = cumsum(filtered);

figure
plot(Vt)
xlabel('Time (ms)')
ylabel('Vt')
hold on

filtered = highpass(data.cage74_2.beforeCre_post.voltage, 150);
Vt = cumsum(filtered);
plot(Vt)
xlabel('Time (ms)')
ylabel('Vt')
hold on

% === Plot inter-breath intervals across time for inspirations and expirations ===