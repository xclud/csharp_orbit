namespace System;

public class PropagationException : Exception
{
    public PropagationException() { }
    public PropagationException(string message) : base(message) { }
    public PropagationException(string message, Exception inner) : base(message, inner) { }
}
