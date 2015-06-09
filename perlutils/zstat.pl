#!/usr/bin/perl
#Zack Ramjan - USC Epigenome Center
#PBS HACK, wrap pbs qstat to look like sge qstat. we emulate array jobs
#by pulling the jobid and tasking from the job name in the form
# name = zXXXX_YYY where XXX is the jobid and YYY is the task id
use Time::Local;
use Term::ANSIColor;
use strict;

my $PBS_QSTAT = "qstat -f";
my $username = $ARGV[0] || $ENV{USER};

#FOR TESTING ONLY
#$PBS_QSTAT = "ssh ramjan\@hpc-login2 qstat -f";
#$username = "ramjan";

my $pbsQstat = `$PBS_QSTAT`;
my @pbsJobs = split("Job Id:", $pbsQstat);

for my $job (@pbsJobs)
{

        my %jobDetails;
        my @jobLines = split("\n", $job);
        chomp $jobLines[0];
	if($jobLines[0] =~ /\s*\d+/)
	{
		$jobDetails{"Job Id"} = $jobLines[0];
		for my $jobProp (@jobLines)
		{
			if($jobProp =~ /\s*(\S+)\s=\s(.+)$/)
			{
				$jobDetails{$1} = $2;
				chomp $jobDetails{$1};
			}
		}


		#clean up pbs output to conform to SGE
		# to simulate arrays we mangle the job names and use it as an id and task
		#JOB ID
		$jobDetails{"Job Id"} =~ /^\s*(\d+)/;
		$jobDetails{"Job Id"} = $1;
		
		my $jobid = my $taskid = "";
		if($jobDetails{"Job_Name"} =~ /z(\d+)_(\d+)/)  
		{
			$jobDetails{"Job Id"} = $1;
			$jobDetails{"Task_Id"} = $2;
		}

		#JOB STATE
		$jobDetails{"job_state"} = lc($jobDetails{"job_state"});
		$jobDetails{"job_state"} =~ s/q/qw/; 

		#JOB START TIME
		my %months =	("Jan"=>1, "Feb"=>2, "Mar"=>3, "Apr"=>4, "May"=>5, "Jun"=>6, 
				"Jul"=>7, "Aug"=>8, "Sep"=>9, "Oct"=>10, "Nov"=>11, "Dec"=>12);
		my @time = split(/\s+/,  $jobDetails{"ctime"}); 
		$jobDetails{"ctime"} = "$months{$time[1]}/$time[2]/$time[4] $time[3]";

		#JOB OWNER
		$jobDetails{"Job_Owner"} =~ s/\@.+$//; 

		#JOB ARRAY

		#JOB NODE
		$jobDetails{"exec_host"} =~ s/(.+?)\/.+/	$1/; 

		#JOB NaME 
		$jobDetails{"Job_Name"} = substr($jobDetails{"Job_Name"},0,15) . "..." .  substr($jobDetails{"Job_Name"},-40) if length($jobDetails{"Job_Name"}) > 59;


		#PRINT RESULTS
		 if ($jobDetails{"Job_Owner"} eq $username || $jobDetails{"queue"} eq $ARGV[0] || $ARGV[0] eq "all")
		{
			print "  " . $jobDetails{"Job Id"} . "\t" ;
			if( $jobDetails{"job_state"} =~ /qw/i) {print color("blue");}
			if( $jobDetails{"job_state"} =~ /r/i) {print color("green");}
			if( $jobDetails{"job_state"} =~ /h/i) {print color("red");}
			if( $jobDetails{"job_state"} =~ /e|c/i) {print color("yellow on_black");}
			print   sprintf("%-60s", $jobDetails{"Job_Name"}) ;
			print color("reset") ;
			print "\t" .  sprintf("%-*s", 10, $jobDetails{"Job_Owner"}) . "\t" . $jobDetails{"job_state"} . "\t";
			print $jobDetails{"ctime"} . "\t" unless $jobDetails{"job_state"} =~ /r/i;
			print " walltime:" . $jobDetails{"resources_used.walltime"} . "\t" if $jobDetails{"job_state"} =~ /r/i;

			print  $jobDetails{"queue"}; 
			if($jobDetails{"exec_host"}) {print color("green");  print  $jobDetails{"exec_host"};}
			else {print color("blue");  print  "\tTBD";}
			print color("reset") . "\t";
			if($jobDetails{"resources_used.mem"}) {print $jobDetails{"resources_used.mem"} . "/";}
			print  $jobDetails{"Resource_List.mem"} . "\n"; 
			print color("reset");

		}


	}
}
