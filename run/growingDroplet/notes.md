# Further test cases and test case variants to test:

## Shrinking droplet
So far, we only tested the growing droplet. However, the same reasonings including
analytical solutions should work for the shrinking droplet. This can test the
receding contact angle.

## Growing droplet with constant contact angle boundary condition
Use the groeing droplet test cases as they are, but with a constant contact
angle boundary condition and compare the interface evolution to the
variant with hysteresis boundary condition.

## Growing droplet with zero inflow
Use zero inflow for the growing droplet with hysteresis boundary condition
and observe the interface evolution. Espiecially, pay attention whether
the shrinking observed initally in case of non-zero inflow continues or
stops in the absence of inflow.

## Growing droplet with original hysteresis and interFoam
Use the original FDT hysteresis boundary condition with interFoam for the
growing droplet cases and check the pinning/spreading behaviour.

Used b01-wetting-benchmark/DropletSpreadingTest/ZeroG as template.

