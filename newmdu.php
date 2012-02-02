<?php

//display errors
ini_set('display_errors',1); 
error_reporting(E_ALL);

//LDAP server setup and variables
$server = "ldap.umd.edu";
$ldap_user = 'cmps-degreeinfo';
$ldap_pass = '8d7F38fWU2ms2vye6Diw9Lri';

//Bind and Search settings
$ldaprdn = 'uid=cmps-degreeinfo,cn=auth,ou=ldap,dc=umd,dc=edu';
//$ldaprdn = "ou=ldap,dc=umd,dc=edu";
$base = "dc=umd, dc=edu";
$filter = "(&(umDepartment=CMNS-Mathematics*)(umEmployee=TRUE)(sn=C*))";
$fields = array("uid","umDisplayName","umNamePrefix","sn","givenName","title","mail","telephoneNumber","postalAddress");

//Connect to the LDAP server
$ldapconn = ldap_connect($server);

//If the connection was successful
if ($ldapconn) {

	//Bind to the LDAP server	
	$ldapbind = ldap_bind($ldapconn,$ldaprdn,$ldap_pass);
	//$ldapbind = ldap_bind($ldapconn,$ldaprdn);
	
	//If the bind was successful
	if ($ldapbind) {
		
		//Search the ldap server
		$search_limit = 0; //set to 0 for all
		$ldapsearch = ldap_search($ldapconn, $base, $filter, $fields,0,$search_limit);
		$ldapSortAttributes = array('givenname','sn');
		foreach ($ldapSortAttributes as $eachSortAttribute) {
			ldap_sort($ldapconn, $ldapsearch, $eachSortAttribute);
		} 		

		//Get the entries
		$ldapinfo = ldap_get_entries($ldapconn, $ldapsearch);
		
		//Result count
		$count = $ldapinfo["count"];
		
		$j=0;

		//If there are entries, put them in an array
		if ($count > 0 ) {
		
			//connect to database
			$mysqli = new mysqli("ub3.ece.umd.edu", "math_joomla", "LLxupAMzYGpUrj7J", "math_joomla");
			if ($mysqli->connect_errno) {
				$error = $mysqli->connect_error;
				exit();
			}
			
			//Delete everyone in category first
			//$mysqli->query("DELETE FROM jos_k2_items WHERE catid = 9 OR catid = 10");
			//Delete all tags
			//$mysqli->query("DELETE FROM jos_k2_tags_xref");
			
			for ($i=0;$i<$count;$i++) {
		
				$itemquery = $mysqli->prepare("SELECT id FROM jos_k2_items WHERE title=? AND trash=0");
				$itemquery->bind_param("s", $name);
				$itemquery->execute();
				if($itemquery->num_rows > 0)
				{
					$itemquery->bind_result($itemID);
					$itemquery->fetch();
					$itemquery->close();
					
					$exists = TRUE;
					
				}
			
				//defaults
				$catid = 8;
				$published = 1;
				$access = 1;
				
									
				//Get ldap results
				$uid = $ldapinfo[$i]["uid"][0];
				
				if (!empty($ldapinfo[$i]["umdisplayname"][0]))
				$name = $ldapinfo[$i]["umdisplayname"][0];
				else $name = "";
				
				if (!empty($ldapinfo[$i]["givenname"][0]))
				$firstname = $ldapinfo[$i]["givenname"][0];
				else $firstname = "";
				
				if (!empty($ldapinfo[$i]["sn"][0]))
				$lastname = $ldapinfo[$i]["sn"][0];
				else $lastname = "";
				
				if (!empty($ldapinfo[$i]["title"][0]))
				$title = $ldapinfo[$i]["title"][0];
				else $title = "";
				
				if (!empty($ldapinfo[$i]["postaladdress"][0]))
				$address = $ldapinfo[$i]["postaladdress"][0];
				else $address = "";
				
				if (!empty($ldapinfo[$i]["telephonenumber"][0]))
				$phone = $ldapinfo[$i]["telephonenumber"][0];
				else $phone = "";
				
				if (!empty($ldapinfo[$i]["mail"][0]))
				$email = $ldapinfo[$i]["mail"][0];
				else $email = "";
				
				if (!empty($ldapinfo[$i]["umnameprefix"][0]))
				$prefix = $ldapinfo[$i]["umnameprefix"][0];
				else $prefix = "";
				
				$out = $name." ".$firstname." ".$lastname." ".$title." ".$address." ".$phone." ".$email;
				echo $out."<br />";
				
				//Get the category
				$firstletter = strtolower(substr($lastname,0,1));
				$result = $mysqli->query("SELECT id FROM jos_k2_categories WHERE alias = '$firstletter'");
				while($row = $result->fetch_assoc()) {
					$catid = $row['id'];
				}
				
				// Prepare the tags
				$tagID = 0;
				$posEmeritus = stripos($title, "Emeritus");
				
				if($posEmeritus !== false)
				{
					echo $name." is emeritus<br />";
					$tagID = 5;
				}
				else
				{
					$posFaculty = stripos($title, "Prof");
					if($posFaculty !== false)
					{
						echo $name." is faculty<br />";
						$tagID = 1;
					}
					$posGrad = stripos($title, "Grad");
					if($posGrad !== false)
					{
						echo $name." is grad<br />";
						$tagID = 2;
					}
					$posLecturer = stripos($title, "Lecturer");
					if($posLecturer !== false)
					{
						echo $name." is lecturer<br />";
						$tagID = 3;
					}
									
				}
				
				// Set up the extra fields
				$extra_fields = "[{\"id\":\"1\",\"value\":\"".$firstname."\"},{\"id\":\"2\",\"value\":\"".$lastname."\"},{\"id\":\"3\",\"value\":\"".$phone."\"},{\"id\":\"4\",\"value\":[\"".$email."\",\"mailto:".$email."\",\"new\"]},{\"id\":\"6\",\"value\":\"".$address."\"},{\"id\":\"7\",\"value\":\"".$title."\"},{\"id\":\"8\",\"value\":\"".$uid."\"},{\"id\":\"9\",\"value\":\"".$prefix."\"},{\"id\":\"10\",\"value\":\"".$name."\"}]";
				$extra_fields_search = $name." ".$firstname." ".$lastname." ".$title." ".$address." ".$phone." ".$email;
				
				//$j++;
				//echo $j." ".$name."<br />";
				
				//insert query
				$query = $mysqli->prepare("INSERT INTO jos_k2_items (title,alias,catid,published,access,extra_fields,extra_fields_search,created) VALUES (?,?,?,?,?,?,?,NOW())");
				$query->bind_param("ssiiiss",$name,$uid,$catid,$published,$access,$extra_fields,$extra_fields_search);
				$query->execute();
				$query->close();
				
				// If we found a tag, insert it
				if($tagID != 0)
				{
					// Now we need to find the query in the database to get the itemID
					// Note that this is a little confusing because in the k2 database, title refers to the full name, while in the LDAP data, it refers to their position at the university.
					$itemquery = $mysqli->prepare("SELECT id FROM jos_k2_items WHERE title=? AND trash=0");
					$itemquery->bind_param("s", $name);
					$itemquery->execute();
					$itemquery->bind_result($itemID);
					$itemquery->fetch();
					$itemquery->close();
				
					echo $itemID." ".$tagID;
					
				
					// We can now insert the tag
					$tagquery = $mysqli->prepare("INSERT INTO jos_k2_tags_xref (tagID,itemID) VALUES (?,?)");
					$tagquery->bind_param("ii",$tagID,$itemID);
					$tagquery->execute();
					printf("Error: %s.\n", $tagquery->error);
					$tagquery->close();
				}
				
								
			}				
		
		}
		
	}
}

?>