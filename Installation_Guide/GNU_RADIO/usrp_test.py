#!/usr/bin/env python3
from gnuradio import gr, uhd, blocks
import time

class usrp_test(gr.top_block):
    def __init__(self):
        gr.top_block.__init__(self, "USRP Test")
        
        # Parameters
        samp_rate = 1e6
        center_freq = 2.4e9
        gain = 30
        
        # Source
        self.usrp_source = uhd.usrp_source(
            ",".join(("", "")),
            uhd.stream_args(cpu_format="fc32", channels=list(range(1))),
        )
        self.usrp_source.set_samp_rate(samp_rate)
        self.usrp_source.set_center_freq(center_freq, 0)
        self.usrp_source.set_gain(gain, 0)
        
        # Sink
        self.blocks_file_sink = blocks.file_sink(
            gr.sizeof_gr_complex*1, 
            '/tmp/usrp_samples.dat', 
            False
        )
        self.connect((self.usrp_source, 0), (self.blocks_file_sink, 0))

if __name__ == '__main__':
    tb = usrp_test()
    tb.start()
    print("Recording for 5 seconds...")
    time.sleep(5)
    tb.stop()
    tb.wait()
    print("Done! Saved to /tmp/usrp_samples.dat")