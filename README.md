# DistributedAberrationCorrection
Distributed Aberration Correction Techniques Based on Tomographic Sound Speed Estimates

Phase aberration correction is the process of applying delay/phase correction to ultrasound signals after applying focusing delays (based on an assumed speed of sound: usually 1540 m/s) to ultrasound channel data prior to summation in a delay-and-sum (DAS) beamformer. Although phase aberration correction is typically caused by sound speed heterogeneities in the medium, the aberration delays/phases are typically modeled independently for each imaging point without reference to the spatial distribution of sound speed in the medium. Here, we show (1) a tomographic sound speed estimator based on phase aberrations after focusing and (2) two different distributed aberration correction techniques that use a tomographic sound speed estimate to create optimally-focused ultrasound images:
* The first distributed aberration correction technique achieves distributed aberration correction by modeling the times-of-flight (TOF) between transducer elements and imaging points by solving the eikonal equation, which relates the sound speed in the medium (based on our tomographic estimate) and the TOF between transducer elements and imaging points. These TOF are used as delays in a DAS beamformer.
* The second distributed aberration correction technique is based on the cross-correlation of transmitted and received wave-fields (see https://github.com/rehmanali1994/FourierDomainBeamformer). This technique models transmitted and receive wave-fields in a heterogeneous sound speed medium by applying the Fourier split-step form of the angular spectrum method. By modeling these wave-fields in a heterogeneous medium, the cross-correlation of transmitted and received wavefields should yield an optimally focused ultrasound image. Note that this imaging approach does not fall under a standard DAS beamforming typically applied in ultrasound imaging.

We provide sample data and algorithms presented in

> Ali, R.; Hyun, D.; Brickson, L.; Dahl, J. "Distributed Aberration Correction Techniques Based on Tomographic Sound Speed Estimates". *Manuscript submitted for publication.*

for the tomographic sound speed estimation technique and the two distributed aberration correction techniques described above. If you use the code/algorithm for research, please cite the above paper.

You can reference a static version of this code by its DOI number: (COMING SOON)

