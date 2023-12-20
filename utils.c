#include <stdio.h>
#include <stdlib.h>

// Warning: if you change the value of MAX_ID_LEN, you might need to change the
// length of the format_specifier string.

#define MAX_ID_LEN 32   // The maximum length of an instance ID
#define MAX_NAME_LEN 64 // The maximum length of an instance name

/*
 * Terminates an AWS instance.
 *
 * @param instance_name The name of the instance to terminate. Max length is 64.
 */
void terminate_aws_instance(const char *instance_name)
{
    FILE *fd;
    // Maximum length of the instance ID (it's actually shorter, but we're being
    // safe for future-proofing)
    char instance_id[MAX_ID_LEN + 1];

    // Length of the AWS command plus the maximum length of the instance name
    char get_id_command[127 + (MAX_NAME_LEN + 1)];

    // Length of the file write command, including the length of the AWS command 
    // and instance name
    char to_file[19 + (127 + MAX_NAME_LEN + 1)];

    // Length of the command to terminate the instance plus the max length of 
    // the instance ID
    char command[58 + (MAX_ID_LEN + 1)];

    // (% + s + null terminator) + length of the instance ID string 
    char format_specifier[3 + 2];

    snprintf(format_specifier, sizeof(format_specifier), "%%%ds", MAX_ID_LEN); // "%ds"

    snprintf(get_id_command, sizeof(get_id_command), "aws ec2 describe-instances --filters \"Name=tag:Name,Values=%s\" --query \"Reservations[*].Instances[*].InstanceId\" --output text", instance_name);

    snprintf(to_file, sizeof(to_file), "%s%s", get_id_command, " > instance_id.txt");

    system(to_file);

    fd = fopen("instance_id.txt", "r");
    if (NULL != fd)
    {
        fscanf(fd, format_specifier, instance_id);

        fclose(fd);

        snprintf(command, sizeof(command), "aws ec2 terminate-instances --instance-ids %s > /dev/null", instance_id);

        system(command);

        remove("instance_id.txt");
    }
}
