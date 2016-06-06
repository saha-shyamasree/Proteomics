use DBI;

my $db="test";
my $server="frontend1:3308";
my $socket="/data/home/btw796/MySQL/thesock2";
my $username="pasa_all";
my $password="pasa_all";

my $dbh=DBI->connect("DBI:mysql:$db:$server:$socket", $username, $password);

unless (ref $dbh) {
        print "Cannot connect to $server: $DBI::errstr";
    }
