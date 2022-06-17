namespace System;

public class PropagationException : Exception
{
    public readonly int Error;
    public PropagationException() { }
    public PropagationException(string message, int error) : base(message)
    {
        Error = error;
    }
    public PropagationException(string message) : base(message) { }
    public PropagationException(string message, Exception inner) : base(message, inner) { }
}
