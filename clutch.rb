#!/usr/bin/ruby

# CLUTCH CONTROL SIMULATION ------------------------------
# by Robin Messenger
#
# This simulates a typical manual transmission vehicle
# and the interaction between a rotating engine, clutch,
# wheels and the road. A genetic algorithm minimizes
# clutch wear.
#
# SIMULATION PARAMETERS
# (All values are in SI units)
$ClutchTolerance=0.02  # clutch slip threshold
$HillGrade=0           # slope of hill
$SpeedGoal=5           # speed goal for GA solutions (m/s)
$MaxSimTime=30         # simulation ended after this (s)
$M=1500                # mass of vehicle (kg)
$I=0.07                # moment of inertia of flywheel
$GearRatio=15.2        # first gear final drive ratio
$WheelRadius=0.24      # radius of wheel (m)
$DT=0.01               # time interval (s)
$MaxTorque=200         # max torque of engine (Nm)
$IdleTorque=136        # engine torque at idle (Nm)
$MaxClutchTorque=400   # max torque before clutch slips (Nm)
$MaxTorqueRPM=2500     # RPM engine prod. max torque
$IdleRPM=800           # RPM at idle
$RedlineRPM=7000       # RPM at redline
$EngineDragTorque=3    # internal engine drag (Nm)
$PopulationSize=100
# --------------------------------------------------------

# These are calculated from above
$MaxTorqueW=$MaxTorqueRPM*Math::PI/60
$IdleW=$IdleRPM*Math::PI/60
$RedlineW=$RedlineRPM*Math::PI/60
$K=$GearRatio/$WheelRadius


# For a given angular velocity, what is engine max torque
def maxTorque w
    if w<$IdleW
        0
    elsif w<$MaxTorqueW
        $IdleTorque+($MaxTorque-$IdleTorque)*
            (w-$IdleW)/
            ($MaxTorqueW-$IdleW)
    else
        $MaxTorque
    end
end


# Represents a simulated driver:
# The gas and clutch pedals are represented as 6th
# degree polynomials
class MyDriver
    attr_accessor :clutch, :gas, :saved_score

    # Initialize random values
    def initialize
        @saved_score=nil
        @clutch=Array[rand-0.5,rand-0.5,rand-0.5,
                      rand-0.5,rand-0.5,rand-0.5]
        @gas=Array[rand-0.5,rand-0.5,rand-0.5,
                   rand-0.5,rand-0.5,rand-0.5]
    end

    # Alter the values randomly
    def mutate
        d=MyDriver.new
        d.clutch=@clutch.map {|c| 
            c*(0.9+0.2*rand)
        }
        d.gas=@gas.map {|c| 
            c*(0.9+0.2*rand)
        }
        d
    end

    # Work out clutch pedal position at time t
    def clutchPedal t
        clutch=0
        @clutch.each_with_index {|c,i| clutch+=c*(t**i)}
        if clutch<0
            0
        elsif clutch<=1
            clutch
        else
            1
        end
    end

    # Work out gas pedal position at time t
    def gasPedal t
        gas=0
        @gas.each_with_index {|c,i| gas+=c*(t**i)}
        if gas<0
            0
        elsif gas<=1
            gas
        else
            1
        end
    end

    # Score driver by running simulation
    def score verbose=false
        return @saved_score if !verbose and @saved_score!=nil

        # Set initial values
        clutchTorque=0
        v=0
        w=$IdleW
        t=0
        energy=0
        a=0
        alpha=0
        force=0
        torque=0
        file=nil
        direction=w-v*$K

        # Open file to write data to
        if verbose
            file=File.open("clutch.dat","w")
            file.puts "#time\tvelocity\tacceleration\t"+
                "RPM x1000\tclutch\tgas"
        end

        # Loop one time step at a time
        while t<$MaxSimTime
            # set angular velocity
            w=v*$K if direction*(w-v*$K)<0  

            direction=w-v*$K
            torque=maxTorque(w)*gasPedal(t)-$EngineDragTorque
            force=-9.8*$HillGrade*$M
            if direction.abs<$ClutchTolerance
                #clutch isn't slipping
                clutchTorque=(torque*$M-$I*force*$K)/($K*$K*$I+$M)
                if clutchTorque.abs>clutchPedal(t)*$MaxClutchTorque
                    if clutchTorque<0
                        clutchTorque=-clutchPedal(t)*
                            $MaxClutchTorque
                    else
                        clutchTorque=clutchPedal(t)*
                            $MaxClutchTorque
                    end
                end
            else
                #clutch is slipping
                clutchTorque=clutchPedal(t)*
                    $MaxClutchTorque
                clutchTorque*=-1 if w<v*$K
            end
            torque-=clutchTorque
            force+=clutchTorque*$K
            alpha=torque/$I
            a=force/$M
            w+=$DT*alpha
            w=0 if w<0
            w=$RedlineW if w>$RedlineW
            v+=$DT*a
            if direction.abs>$ClutchTolerance
                energy+=clutchPedal(t)*
                        $MaxClutchTorque*
                        (w-v*$K).abs*
                        $DT
            end
            file.puts "#{t}\t"+  # Print parameters to a file
                      "#{v}\t"+  # if in verbose mode
                      "#{a}\t"+
                      "#{w*60/Math::PI/1000}\t"+
                      "#{clutchPedal(t)}\t"+
                      "#{gasPedal(t)}" if verbose
            t+=$DT
            if w<$IdleW
                if verbose
                    puts "Car stalled, maximum speed:"+
                        " #{v} meters/sec"
                    file.close 
                end
                @saved_score=-1000000+v
                return @saved_score
            end
            if v>=$SpeedGoal
                if verbose
                    puts "Started successfully, "+
                        "clutch dissipated #{energy} J"
                    file.close 
                end
                @saved_score=1/(1+energy)
                return @saved_score
            end
        end
        if verbose
            puts "Failed to start in time, "+
                "#{(w-v*$K).abs*60/Math::PI} RPM mismatch"
            file.close
        end
        @saved_score=-(w-v*$K).abs
        return @saved_score
    end
end


# Create a population
drivers=Array.new
$PopulationSize.times {drivers.push MyDriver.new}
best=MyDriver.new

loop do
    $PopulationSize.times do
        #pick two individuals
        a=rand(drivers.size)
        b=(a+1)%drivers.size

        #replace worst with best
        if drivers[a].score>drivers[b].score
            drivers[b]=drivers[a].mutate
            best=drivers[a] if drivers[a].score>best.score
        else
            drivers[a]=drivers[b].mutate
            best=drivers[b] if drivers[b].score>best.score
        end
    end

    #save data from best one
    best.score(true)
end

