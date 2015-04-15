/**
 * A very simple options parser.  Only allows one name per option, returns everything as strings.
 * Takes a string[] of option names, and a string[] to return the values for the options.  Options that require
 * no parameters are returned as the string "true"
 * 
 * 
 * @author matthewb
 *
 */

public class ParseOpts
{

	/**
	 * @param args
	 */
	public static void main(String[] args)
	{
		String[] on = {"A", "S", "D", "F"};
		String[] vals = new String[on.length];
		String warns = ParseOpts.parse(on, args, vals);
		if(warns.length() > 0)
		{
			System.out.println("Warning: ");
			System.out.println(warns);
		}
		for(int i = 0; i < on.length; i++)
		{
			System.out.println(on[i]+" "+vals[i]);
		}

	}
	
	/**
	 * Parses args based on optnames, into values.
	 * @param optnames The names of allowable options
	 * @param args The input arguments
	 * @param values The return values for each optname (co-indexed)
	 * @return
	 */
	public static String parse(String[] optnames, String[] args, String[] values)
	{
		StringBuffer warnings = new StringBuffer(100);
		for(int i = 0; i < args.length; i++)
		{
			String curr = args[i];
			if(curr.charAt(0) == '-')
			{
				curr = curr.substring(1);
				String value = null;
				String optname = null;
				if(curr.indexOf('=')!= -1)
				{
					optname = curr.substring(0, curr.indexOf('='));
					value = curr.substring(curr.indexOf('=')+1);
				}
				else
				{
					if(i<args.length -1)
					{
						String next = args[i+1];
						if(next.charAt(0) == '-')
						{
							optname = curr;
							value = "true";
						}
						else
						{
							optname = curr;
							value = next;
							i++;
						}	
					}
					else
					{
						optname = curr;
						value = "true";						
					}
				}
				boolean found = false;
				
				for(int j = 0; j < optnames.length; j++)
				{
					if(optname.equals(optnames[j]))
					{
						found = true;
						values[j] = value;
						break;
					}	
				}
				if(!found)
				{
					warnings.append("No such option: "+optname+"\n");
				}
			}
			else
			{
				warnings.append("Badly formed option: "+curr+"\n");
			}
		}
		return warnings.toString();
	}
	

}
