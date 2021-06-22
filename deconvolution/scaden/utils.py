
import sys
import logging
import subprocess

logger = logging.getLogger(__name__)


def _check_call_step(cmd, current_step = -1, init_step = -1, call=True, 
        raise_on_error=True):
    
    logging.info(cmd)
    ret_code = 0

    if current_step >= init_step:
        if call:
            #logging.info(cmd)
            logging.info("calling")
            ret_code = subprocess.call(cmd, shell=True)

            if raise_on_error and (ret_code != 0):
                raise subprocess.CalledProcessError(ret_code, cmd)
            elif (ret_code != 0):
                msg = ("The command returned a non-zero return code\n\t{}\n\t"
                    "Return code: {}".format(cmd, ret_code))
                logger.warning(msg)
        else:
            msg = "skipping due to --do-not-call flag"
            logging.info(msg)
    else:
        msg = "skipping due to --init-step; {}, {}".format(current_step, init_step)
        logging.info(msg)

    return ret_code


def _check_call(cmd, call=True, raise_on_error=True):
    return _check_call_step(cmd, call=call, raise_on_error=raise_on_error)

def _check_output_step(cmd, current_step = 0, init_step = 0, raise_on_error=True):

    logging.info(cmd)
    if current_step >= init_step:
        logging.info("calling")

        try:
            out = subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as exc:
            if raise_on_error:
                raise exc
        return out.decode()

def _check_output(cmd, call=True, raise_on_error=True):
    current_step = 1
    init_step = 1
    if not call:
        init_step = 2
    return _check_output_step(cmd, current_step, init_step, 
        raise_on_error=raise_on_error)


def _check_srun(cmd, call=True):
    cmd = 'srun {}'.format(cmd)
    _check_call(cmd, call)

def _check_sbatch(cmd, call=True, num_cpus=1, mem="25G", time=None, 
        partitions=None, dependencies=None, no_output=False, no_error=False, 
        use_slurm=False, mail_type=['FAIL', 'TIME_LIMIT'], mail_user=None,
        stdout_file=None, stderr_file=None,
        args=None, export=None, num_gpus=None, gpu_name=None):

    # use args if they are present
    # export (dict) not in args
    if args is not None:
        call = not args.do_not_call
        num_cpus = args.num_cpus
        mem = args.mem
        time = args.time
        partitions = args.partitions
        no_output = args.no_output
        no_error = args.no_error
        use_slurm = args.use_slurm
        mail_type = args.mail_type
        mail_user = args.mail_user
        stdout_file = args.stdout_file
        stderr_file = args.stderr_file

    output_str = "--output=slurm-%J.out"
    if stdout_file is not None:
        output_str = "--output={}".format(stdout_file)
    if no_output:
        output_str = "--output=/dev/null"

    error_str = "--error=slurm-%J.err" 
    if stderr_file is not None:
        error_str = "--error={}".format(stderr_file)
    if no_error:
        error_str = "--error=/dev/null"

    dependencies_str = ""
    if dependencies is not None:
        dependencies_str = ':'.join(str(d) for d in dependencies)
        dependencies_str = "--dependency=afterok:{}".format(dependencies_str)
    
    export_str = ""
    if export is not None:
        export_str = ','.join(str(d) for d in [f"{k}={v}" for k,v in export.items()])
        export_str = "--export=ALL,{}".format(export_str)
        
    # check if we actually want to use SLURM
    msg = "_check_sbatch.use_slurm: {}, call: {}".format(use_slurm, call)
    logger.debug(msg)

    # anyway, make sure to remove the --use-slurm option
    cmd = cmd.replace("--use-slurm", "")
    
    gres_str = ""
    if num_gpus is not None and gpu_name is not None:
        gres_str = "--gres=gpu:{}:{}".format(gpu_name, num_gpus)
    
    if use_slurm:
        time_str = ""
        if time is not None:
            time_str = "--time {}".format(time)

        mem_str = ""
        if mem is not None:
            mem_str = "--mem={}".format(mem)

        partitions_str = ""
        if partitions is not None:
            partitions_str = "-p {}".format(partitions)

        num_cpus_str = ""
        if num_cpus is not None:
            num_cpus_str = "--cpus-per-task {}".format(num_cpus)

        mail_type_str = ""
        if mail_type is not None:
            mail_type_str = "--mail-type {}".format(','.join(mail_type))

        mail_user_str = ""
        if mail_user is not None:
            mail_user_str = "--mail-user {}".format(mail_user)
        else:
            # if we did not give a mail user, then do not specify the mail types
            mail_type_str = ""

        cmd = ("sbatch {} {} --ntasks 1 {} {} "
            "{} {} {} {} {} {} {} {}".format(time_str, mem_str, partitions_str, num_cpus_str, gres_str,
                                             dependencies_str, output_str, error_str, mail_type_str, 
                                             mail_user_str, export_str, cmd))

        output = _check_output(cmd, call=call)

        # and parse out the job id
        if call:
            job_id = output.strip().split()[-1]
        else:
            job_id = None
        return job_id
    else:
        _check_call(cmd, call=call)
        return None

mail_type_choices = [
    'NONE', 
    'BEGIN', 
    'END', 
    'FAIL', 
    'REQUEUE', 
    'ALL',
    'STAGE_OUT',
    'TIME_LIMIT', 
    'TIME_LIMIT_90',
    'TIME_LIMIT_80',
    'TIME_LIMIT_50',
    'ARRAY_TASKS'
]

def _add_sbatch_options(parser, num_cpus=1, mem="2G", time=None, 
            stdout_file=None, stderr_file=None,
            partitions=None, mail_type=['FAIL', 'TIME_LIMIT'], mail_user=None):

    slurm_options = parser.add_argument_group("slurm options")

    slurm_options.add_argument('--num-cpus', help="The number of CPUs to use", type=int, 
        default=num_cpus)
    slurm_options.add_argument('--mem', help="The amount of RAM to request", default=mem)
    slurm_options.add_argument('--time', help="The amount of time to request", default=time)
    slurm_options.add_argument('--partitions', help="The partitions to request", default=partitions)
    slurm_options.add_argument('--no-output', help="If this flag is present, stdout "
        "will be redirected to /dev/null", action='store_true')
    slurm_options.add_argument('--no-error', help="If this flag is present, stderr "
        "will be redirected to /dev/null", action='store_true')
    slurm_options.add_argument('--stdout-file', help="If this is present and the "
        "--no-output flag is not given, then stdout will be directed to this "
        "file rather than slurm-<job>.out", default=stdout_file)
    slurm_options.add_argument('--stderr-file', help="If this is present and the "
        "--no-error flag is not given, then stderr will be directed to this "
        "file rather than slurm-<job>.err", default=stderr_file)
    slurm_options.add_argument('--do-not-call', help="If this flag is present, then the commands "
        "will not be executed (but will be printed).", action='store_true')
    slurm_options.add_argument('--use-slurm', help="If this flag is present, the program calls "
        "will be submitted to SLURM.", action='store_true')
    slurm_options.add_argument('--mail-type', help="When to send an email notifcation "
        "of the job status. See official documentation for a description of the "
        "values. If a mail-user is not specified, this will revert to 'None'.", 
        nargs='*', choices=mail_type_choices, default=mail_type)
    slurm_options.add_argument('--mail-user', help="To whom an email will be sent, in "
        "accordance with mail-type", default=mail_user)

def _get_slurm_options_string(args):

    args_dict = vars(args)

    # first, pull out the text arguments
    sbatch_options = ['num_cpus', 'mem', 'time', 'partitions', 'mail_user', 
        'stdout_file', 'stderr_file']

    # create a new dictionary mapping from the flag to the value
    sbatch_flags_and_vals = {'--{}'.format(o.replace('_', '-')) : args_dict[o] 
        for o in sbatch_options if args_dict[o] is not None}

    # and we need to handle mail-type
    mail_types = " ".join(mt for mt in args_dict['mail_type'])
    mail_types = "--mail-type {}".format(mail_types)

    s = ' '.join("{} '{}'".format(k,v) for k,v in sbatch_flags_and_vals.items())

    s = "{} {}".format(mail_types, s)

    # and check the flags
    if args.no_output:
        s = "--no-output {}".format(s)

    if args.no_error:
        s = "--no-error {}".format(s)

    if args.do_not_call:
        s = "--do-not-call {}".format(s)

    if args.use_slurm:
        s = "--use-slurm {}".format(s)

    return s


def _add_logging_options(parser, default_log_file=""):

    logging_options = parser.add_argument_group("logging options")

    default_log_file = ""
    logging_level_choices = ['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
    default_logging_level = 'WARNING'
    default_specific_logging_level = 'NOTSET'

    logging_options.add_argument('--log-file', help="This option specifies a file to "
        "which logging statements will be written (in addition to stdout and "
        "stderr, if specified)", default=default_log_file)
    logging_options.add_argument('--log-stdout', help="If this flag is present, then "
        "logging statements will be written to stdout (in addition to a file "
        "and stderr, if specified)", action='store_true')
    logging_options.add_argument('--no-log-stderr', help="Unless this flag is present, then "
        "logging statements will be written to stderr (in addition to a file "
        "and stdout, if specified)", action='store_true')

    logging_options.add_argument('--logging-level', help="If this value is specified, "
        "then it will be used for all logs", choices=logging_level_choices,
        default=default_logging_level)
    logging_options.add_argument('--file-logging-level', help="The logging level to be "
        "used for the log file, if specified. This option overrides "
        "--logging-level.", choices=logging_level_choices, 
        default=default_specific_logging_level)
    logging_options.add_argument('--stdout-logging-level', help="The logging level to be "
        "used for the stdout log, if specified. This option overrides "
        "--logging-level.", choices=logging_level_choices, 
        default=default_specific_logging_level)
    logging_options.add_argument('--stderr-logging-level', help="The logging level to be "
        "used for the stderr log, if specified. This option overrides "
        "--logging-level.", choices=logging_level_choices, 
        default=default_specific_logging_level)
    
    
def _update_logging(args, logger=None, 
        format_str='%(levelname)-8s %(name)-8s %(asctime)s : %(message)s'):

    # find the root logger if another logger is not specified
    if logger is None:
        logger = logging.getLogger('')
            
    logger.handlers = []

    # set the base logging level
    level = logging.getLevelName(args.logging_level)
    logger.setLevel(level)

    # now, check the specific loggers

    if len(args.log_file) > 0:
        h = logging.FileHandler(args.log_file)
        formatter = logging.Formatter(format_str)
        h.setFormatter(formatter)
        if args.file_logging_level != 'NOTSET':
            l = logging.getLevelName(args.file_logging_level)
            h.setLevel(l)
        logger.addHandler(h)

    if args.log_stdout:
        h = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter(format_str)
        h.setFormatter(formatter)
        if args.stdout_logging_level != 'NOTSET':
            l = logging.getLevelName(args.stdout_logging_level)
            h.setLevel(l)
        logger.addHandler(h)

    log_stderr = not args.no_log_stderr
    if log_stderr:
        h = logging.StreamHandler(sys.stderr)
        formatter = logging.Formatter(format_str)
        h.setFormatter(formatter)
        if args.stderr_logging_level != 'NOTSET':
            l = logging.getLevelName(args.stderr_logging_level)
            h.setLevel(l)
        logger.addHandler(h)


def _check_programs_exist(programs, raise_on_error=True, package_name=None, 
            logger=logger):

    import shutil

    missing_programs = []
    for program in programs:
        exe_path = shutil.which(program)

        if exe_path is None:
            missing_programs.append(program)

    if len(missing_programs) > 0:
        missing_programs = ' '.join(missing_programs)
        msg = "The following programs were not found: " + missing_programs

        if package_name is not None:
            msg = msg + ("\nPlease ensure the {} package is installed."
                .format(package_name))

        if raise_on_error:
            raise EnvironmentError(msg)
        else:
            logger.warning(msg)

    return missing_programs

