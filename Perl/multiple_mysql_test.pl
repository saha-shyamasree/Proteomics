use Carp;
use DBI;

my $db="G79";
my $server="frontend1";
my $socket="mysql_socket=/tmp/mysql.sock";
my $username="pasa_readonly";
my $password="pasa_readonly";
my $dbh = DBI->connect("DBI:mysql:$db:$server:$socket", $username, $password);
unless (ref $dbh) {
	croak "Cannot connect to $server: $DBI::errstr";
}
$dbh->{RaiseError} = 1; #tu
